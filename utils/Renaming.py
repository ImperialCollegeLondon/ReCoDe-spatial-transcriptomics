#!/usr/bin/env python3
# Usage:
# python3 Renaming.py -d /path/to/the/folder/where/your/sample/folders/are -m path/to/the/mapping/file --rename

# You need to create a Mapping.md file, like this:
# ch0000_dapi.ome.tif: morphology_focus_0000.ome.tif
# ch0001_atp1a1_cd45_e-cadherin.ome.tif: morphology_focus_0001.ome.tif
# ch0002_18s.ome.tif: morphology_focus_0002.ome.tif
# ch0003_alphasma_vimentin.ome.tif: morphology_focus_0003.ome.tif

"""Rename morphology_focus OME-TIFFs and update OME-XML in TIFF comments.
This script reads a simple mapping file (default: data/xenium_output_format.md)
containing lines like:
  morphology_focus_0000.ome.tif: DAPI image
It will scan the `morphology_focus` folders under the data directory, plan
renames based on the mapping, optionally perform backups, rename files, and
update the OME-XML stored in the TIFF comment to (1) set a sensible Channel
Name and (2) update any occurrences of the old filename to the new filename
inside the XML.
This is conservative: it only renames files that it finds mappings for and
avoids clobbering targets unless `--overwrite` is used.
"""

from pathlib import Path
import argparse
import logging
import re
import os
import shutil
import sys
from xml.etree import ElementTree as ET


def parse_mapping_file(path: Path):
    mapping = {}
    with path.open("r", encoding="utf8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("```"):
                continue
            if ":" not in line:
                continue
            key, val = line.split(":", 1)
            mapping[key.strip()] = val.strip()
    return mapping


def find_morphology_tifs(data_dir: Path):
    # look for any morphology_focus directories and collect *.ome.tif files
    files = []
    for p in data_dir.rglob("morphology_focus"):
        if p.is_dir():
            files.extend(sorted([f for f in p.glob("*.ome.tif") if f.is_file()]))
    return files


def slugify(text: str) -> str:
    text = text.lower()
    text = re.sub(r"[\s/]+", "_", text)
    text = re.sub(r"[^a-z0-9_]+", "", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text or "unnamed"


def build_rename_map(files, mapping):
    rename_map = {}
    used = set()
    for f in files:
        desc = mapping.get(f.name) or mapping.get(f.stem)
        print(desc)
        if not desc:
            # try to extract an index from patterns like ch0000 or morphology_focus_0000
            m = re.search(r"(?:ch|morphology_focus)_0*(\d+)", f.name, flags=re.IGNORECASE)
            if m:
                idx = int(m.group(1))
                key = f"morphology_focus_{idx:04d}.ome.tif"
                desc = mapping.get(key)

        if not desc:
            continue

        # If the mapping value is already a filename (e.g. morphology_focus_0000.ome.tif),
        # use it directly. Otherwise, slugify the description and preserve suffix.
        if re.search(r"\.ome\.tif$", desc, flags=re.IGNORECASE):
            # ensure we only take the basename portion
            candidate = Path(desc).name
            base_part = candidate.rsplit('.', 2)[0] if '.' in candidate else candidate
            suffix = ''.join(f.suffixes) if f.suffixes else f.suffix
            # if candidate already ends with the same suffix (like .ome.tif), use it
            if candidate.lower().endswith(suffix.lower()):
                new_name = candidate
            else:
                new_name = f"{candidate}"
        else:
            base = slugify(desc)
            suffix = ''.join(f.suffixes) if f.suffixes else f.suffix
            new_name = f"{base}{suffix}"
        i = 1
        while new_name in used or (f.with_name(new_name)).exists():
            new_name = f"{base}_{i}{suffix}"
            i += 1

        used.add(new_name)
        rename_map[f] = f.with_name(new_name)
    return rename_map


def extract_label_from_key(key: str) -> str:
    """Derive a human-friendly label from a mapping key (e.g. ch0000_dapi.ome.tif -> 'DAPI')."""
    name = Path(key).stem
    m = re.match(r"ch0*(\d+)_(.+)", name, flags=re.IGNORECASE)
    if m:
        label = m.group(2)
    else:
        # fallback: remove numeric prefix like morphology_focus_0000
        m2 = re.match(r"(?:morphology_focus|ch)[_0-9]*(?:_)?(.+)", name, flags=re.IGNORECASE)
        label = m2.group(1) if m2 else name
    # replace underscores with spaces and title-case
    label = label.replace('_', ' ').strip()
    # keep uppercase acronyms if present by heuristic (all letters uppercase)
    if label.isupper():
        return label
    return label.title()


def set_channel_name_and_replace_xml(xml_text: str, channel_name: str, replacements: dict) -> str:
    # set Channel Name (conservative: first Channel element)
    try:
        root = ET.fromstring(xml_text)
    except Exception:
        # if parsing fails, still attempt textual replacements
        out = xml_text
        for old, new in replacements.items():
            out = out.replace(old, new)
        return out

    channels = [el for el in root.iter() if el.tag.split('}')[-1].lower() == 'channel']
    if channels:
        ch = channels[0]
        ch.set('Name', channel_name)
        # set or create Name child
        name_child = None
        for child in ch:
            if child.tag.split('}')[-1].lower() == 'name':
                name_child = child
                break
        if name_child is None:
            # preserve namespace if any
            if '}' in ch.tag:
                ns = ch.tag.split('}')[0].strip('{')
                name_tag = f"{{{ns}}}Name"
            else:
                name_tag = 'Name'
            name_child = ET.Element(name_tag)
            name_child.text = channel_name
            ch.insert(0, name_child)
        else:
            name_child.text = channel_name

    new_xml = ET.tostring(root, encoding='utf-8')
    if isinstance(new_xml, bytes):
        new_xml = new_xml.decode('utf-8')

    # Walk the entire tree and replace filename occurrences in attributes and text
    for el in root.iter():
        # attributes
        for attr, val in list(el.attrib.items()):
            new_val = val
            for old, new in replacements.items():
                if old in new_val:
                    new_val = new_val.replace(old, new)
            if new_val != val:
                el.set(attr, new_val)

        # text
        if el.text:
            new_text = el.text
            for old, new in replacements.items():
                if old in new_text:
                    new_text = new_text.replace(old, new)
            if new_text != el.text:
                el.text = new_text

        # tail
        if el.tail:
            new_tail = el.tail
            for old, new in replacements.items():
                if old in new_tail:
                    new_tail = new_tail.replace(old, new)
            if new_tail != el.tail:
                el.tail = new_tail

    # final serialization
    new_xml = ET.tostring(root, encoding='utf-8')
    if isinstance(new_xml, bytes):
        new_xml = new_xml.decode('utf-8')

    return new_xml


def write_ascii_tiff_comment(path: Path, xml_text: str):
    from tifffile import tiffcomment

    try:
        ascii_xml = xml_text.encode('ascii', 'xmlcharrefreplace').decode('ascii')
    except Exception:
        ascii_xml = xml_text

    tiffcomment(str(path), comment=ascii_xml)


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument('-d', '--data-dir', default='./data')
    p.add_argument('-m', '--mapping', default='./data/xenium_output_format.md')
    p.add_argument('--rename', action='store_true')
    p.add_argument('--backup', action='store_true')
    p.add_argument('--dry-run', action='store_true')
    p.add_argument('--overwrite', action='store_true')
    p.add_argument('-v', '--verbose', action='store_true')
    args = p.parse_args(argv)

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format='%(levelname)s: %(message)s')

    for subdir_name in os.listdir(os.path.abspath(args.data_dir)):
        subdir = os.path.join(os.path.abspath(args.data_dir), subdir_name)

        data_dir = Path(os.path.abspath(subdir))
        if not data_dir.exists():
            logging.error('Data directory does not exist: %s', data_dir)
            sys.exit(2)

        mapping_file = Path(args.mapping)
        if not mapping_file.exists():
            logging.error('Mapping file not found: %s', mapping_file)
            sys.exit(2)

        mapping = parse_mapping_file(mapping_file)
        print(mapping)
        files = find_morphology_tifs(data_dir)
        if not files:
            logging.info('No morphology_focus *.ome.tif files found under %s', data_dir)
            return

        logging.info('Found %d morphology_focus *.ome.tif files', len(files))

        rename_map = {}
        if args.rename:
            rename_map = build_rename_map(files, mapping)
            print(rename_map)
            if not rename_map:
                logging.warning('No files matched mapping enstries; nothing to rename')

            else:
                logging.info('Planned renames:')
                for old, new in rename_map.items():
                    logging.info('  %s -> %s', old.name, new.name)

                if args.dry_run:
                    logging.info('Dry-run: not performing renames')
                else:
                    for old, new in rename_map.items():
                        if new.exists() and not args.overwrite:
                            logging.warning('Target exists and --overwrite not set: %s -> %s', old, new)
                            continue
                        if args.backup:
                            bak = old.with_suffix(old.suffix + '.bak')
                            shutil.copy2(old, bak)
                            logging.info('Backed up %s -> %s', old.name, bak.name)
                        logging.info('Renaming %s -> %s', old.name, new.name)
                        old.replace(new)

                    # refresh file list
                    files = find_morphology_tifs(data_dir)

        # build replacements map: include mapping file (old basename -> new basename)
        replacements = {}
        for k, v in mapping.items():
            replacements[Path(k).name] = Path(v).name
        # include any renames performed earlier
        for old, new in rename_map.items():
            replacements[old.name] = new.name

        # inverse map from mapping values (targets) to mapping keys (sources)
        mapping_value_to_key = {Path(v).name: Path(k).name for k, v in mapping.items()}

        # process each file: set channel name and update xml references
        successes = 0
        for f in files:
            # determine desired channel name: find mapping key (source) for this target file
            key = mapping_value_to_key.get(f.name)
            desc = None
            if key:
                # derive a readable label from the mapping key
                label = extract_label_from_key(key)
                desc = label

            # fallback: if mapping provided a direct mapping for this filename
            if desc is None:
                desc = mapping.get(f.name) or mapping.get(f.stem)

            if not desc:
                logging.info('No mapping for %s; skipping XML edits', f.name)
                continue

            if args.dry_run:
                logging.info('Dry-run: would set channel name for %s -> %s', f.name, desc)
                successes += 1
                continue

            from tifffile import tiffcomment

            try:
                comment = tiffcomment(str(f))
            except Exception as exc:
                logging.error('Failed to read TIFF comment from %s: %s', f, exc)
                continue

            if not comment:
                logging.warning('No TIFF comment in %s; skipping', f)
                continue

            xml_text = comment.decode('utf-8', errors='replace') if isinstance(comment, bytes) else str(comment)
            new_xml = set_channel_name_and_replace_xml(xml_text, desc, replacements)

            try:
                # backup original TIFF comment by copying file if requested
                if args.backup:
                    bak = f.with_suffix(f.suffix + '.pre_xml_bak')
                    shutil.copy2(f, bak)
                    logging.info('Backed up %s -> %s', f.name, bak.name)

                write_ascii_tiff_comment(f, new_xml)
                logging.info('Updated XML for %s', f.name)
                successes += 1
            except Exception as exc:
                logging.exception('Failed to write TIFF comment for %s: %s', f, exc)

        logging.info('Completed: %d/%d files updated', successes, len(files))


if __name__ == '__main__':
    main()
