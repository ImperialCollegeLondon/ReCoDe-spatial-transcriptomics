#!/usr/bin/env python3
"""Batch rename morphology_focus OME-TIFFs across multiple Xenium output directories.

Usage:
    python3 xenium_batch_rename.py \
        -d /path/to/xenium/run/directory \
        -m path/to/mapping/file --rename

This script will:
1. Find all output-* directories in the specified data directory
2. For each output directory, process the morphology_focus folder
3. Rename files according to the mapping file
4. Update OME-XML in TIFF comments with proper channel names

You need to create a mapping file like xenium_mapping.md:
    ch0000_dapi.ome.tif: morphology_focus_0000.ome.tif
    ch0001_atp1a1_cd45_e-cadherin.ome.tif: morphology_focus_0001.ome.tif
    ch0002_18s.ome.tif: morphology_focus_0002.ome.tif
    ch0003_alphasma_vimentin.ome.tif: morphology_focus_0003.ome.tif
"""

import argparse
import logging
import re
import shutil
import sys
from pathlib import Path
from xml.etree import ElementTree as ET


def parse_mapping_file(path: Path):
    """Parse the mapping file and return a dictionary."""
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
    """Find all *.ome.tif files in morphology_focus directories."""
    files = []
    for p in data_dir.rglob("morphology_focus"):
        if p.is_dir():
            files.extend(sorted([f for f in p.glob("*.ome.tif") if f.is_file()]))
    return files


def slugify(text: str) -> str:
    """Convert text to a filesystem-safe slug."""
    text = text.lower()
    text = re.sub(r"[\s/]+", "_", text)
    text = re.sub(r"[^a-z0-9_]+", "", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text or "unnamed"


def build_rename_map(files, mapping):
    """Build a mapping of old filenames to new filenames."""
    rename_map = {}
    used = set()

    for f in files:
        desc = mapping.get(f.name) or mapping.get(f.stem)

        if not desc:
            # try to extract an index from patterns like ch0000 or morphology_focus_0000
            m = re.search(
                r"(?:ch|morphology_focus)_0*(\d+)", f.name, flags=re.IGNORECASE
            )
            if m:
                idx = int(m.group(1))
                key = f"morphology_focus_{idx:04d}.ome.tif"
                desc = mapping.get(key)

        if not desc:
            continue

        suffix = "".join(f.suffixes) if f.suffixes else f.suffix

        # If the mapping value is already a filename (morphology_focus_0000.ome.tif)
        # use it directly. Otherwise, slugify the description and preserve suffix.
        if re.search(r"\.ome\.tif$", desc, flags=re.IGNORECASE):
            candidate = Path(desc).name
            base = Path(candidate).stem  # define base even for .ome.tif cases
            if candidate.lower().endswith(suffix.lower()):
                new_name = candidate
            else:
                new_name = f"{candidate}"
        else:
            base = slugify(desc)
            new_name = f"{base}{suffix}"

        # Handle naming conflicts
        i = 1
        while new_name in used or (f.with_name(new_name)).exists():
            new_name = f"{base}_{i}{suffix}"
            i += 1

        used.add(new_name)
        rename_map[f] = f.with_name(new_name)

    return rename_map


def extract_label_from_key(key: str) -> str:
    """Derive a human-friendly label from a mapping key."""
    name = Path(key).stem
    m = re.match(r"ch0*(\d+)_(.+)", name, flags=re.IGNORECASE)
    if m:
        label = m.group(2)
    else:
        # fallback: remove numeric prefix like morphology_focus_0000
        m2 = re.match(
            r"(?:morphology_focus|ch)[_0-9]*(?:_)?(.+)", name, flags=re.IGNORECASE
        )
        label = m2.group(1) if m2 else name

    # replace underscores with spaces and title-case
    label = label.replace("_", " ").strip()
    if label.isupper():
        return label
    return label.title()


def set_channel_name_and_replace_xml(
    xml_text: str, channel_name: str, replacements: dict
) -> str:
    """Set channel name and replace filename references in OME-XML."""
    try:
        root = ET.fromstring(xml_text)
    except Exception:
        # if parsing fails, still attempt textual replacements
        out = xml_text
        for old, new in replacements.items():
            out = out.replace(old, new)
        return out

    # Set Channel Name (conservative: first Channel element)
    channels = [el for el in root.iter() if el.tag.split("}")[-1].lower() == "channel"]
    if channels:
        ch = channels[0]
        ch.set("Name", channel_name)
        # set or create Name child
        name_child = None
        for child in ch:
            if child.tag.split("}")[-1].lower() == "name":
                name_child = child
                break
        if name_child is None:
            # preserve namespace if any
            if "}" in ch.tag:
                ns = ch.tag.split("}")[0].strip("{")
                name_tag = f"{{{ns}}}Name"
            else:
                name_tag = "Name"
            name_child = ET.Element(name_tag)
            name_child.text = channel_name
            ch.insert(0, name_child)
        else:
            name_child.text = channel_name

    # Walk the entire tree and replace filename occurrences
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
    new_xml = ET.tostring(root, encoding="utf-8")
    if isinstance(new_xml, bytes):
        new_xml = new_xml.decode("utf-8")

    return new_xml


def write_ascii_tiff_comment(path: Path, xml_text: str):
    """Write ASCII XML to TIFF comment."""
    from tifffile import tiffcomment

    try:
        ascii_xml = xml_text.encode("ascii", "xmlcharrefreplace").decode("ascii")
    except Exception:
        ascii_xml = xml_text

    tiffcomment(str(path), comment=ascii_xml)


def process_output_directory(output_dir: Path, mapping: dict, args):
    """Process a single output directory."""
    logging.info(f"Processing directory: {output_dir.name}")

    if not output_dir.exists():
        logging.error(f"Output directory does not exist: {output_dir}")
        return False

    files = find_morphology_tifs(output_dir)
    if not files:
        logging.info(f"No morphology_focus *.ome.tif files found in {output_dir}")
        return True

    logging.info(
        f"Found {len(files)} morphology_focus *.ome.tif files in {output_dir.name}"
    )

    rename_map = {}
    if args.rename:
        rename_map = build_rename_map(files, mapping)
        if not rename_map:
            logging.warning(
                f"No files matched mapping entries in {output_dir.name};"
                f"nothing to rename"
            )
        else:
            logging.info(f"Planned renames for {output_dir.name}:")
            for old, new in rename_map.items():
                logging.info("  %s -> %s", old.name, new.name)

            if not args.dry_run:
                for old, new in rename_map.items():
                    if new.exists() and not args.overwrite:
                        logging.warning(
                            "Target exists and --overwrite not set: %s -> %s", old, new
                        )
                        continue
                    if args.backup:
                        bak = old.with_suffix(old.suffix + ".bak")
                        shutil.copy2(old, bak)
                        logging.info("Backed up %s -> %s", old.name, bak.name)
                    logging.info("Renaming %s -> %s", old.name, new.name)
                    old.replace(new)

                # refresh file list after renaming
                files = find_morphology_tifs(output_dir)

    # Build replacements map for XML updates
    replacements = {}
    for k, v in mapping.items():
        replacements[Path(k).name] = Path(v).name
    for old, new in rename_map.items():
        replacements[old.name] = new.name

    # Inverse map from mapping values to mapping keys
    mapping_value_to_key = {Path(v).name: Path(k).name for k, v in mapping.items()}

    # Process each file: set channel name and update xml references
    successes = 0
    for f in files:
        # Determine desired channel name
        key = mapping_value_to_key.get(f.name)
        desc = None
        if key:
            label = extract_label_from_key(key)
            desc = label

        # Fallback: if mapping provided a direct mapping for this filename
        if desc is None:
            desc = mapping.get(f.name) or mapping.get(f.stem)

        if not desc:
            logging.info("No mapping for %s; skipping XML edits", f.name)
            continue

        if args.dry_run:
            logging.info("Dry-run: would set channel name for %s -> %s", f.name, desc)
            successes += 1
            continue

        from tifffile import tiffcomment

        try:
            comment = tiffcomment(str(f))
        except Exception as exc:
            logging.error("Failed to read TIFF comment from %s: %s", f, exc)
            continue

        if not comment:
            logging.warning("No TIFF comment in %s; skipping", f)
            continue

        xml_text = (
            comment.decode("utf-8", errors="replace")
            if isinstance(comment, bytes)
            else str(comment)
        )
        new_xml = set_channel_name_and_replace_xml(xml_text, desc, replacements)

        try:
            # backup original TIFF comment by copying file if requested
            if args.backup:
                bak = f.with_suffix(f.suffix + ".pre_xml_bak")
                shutil.copy2(f, bak)
                logging.info("Backed up %s -> %s", f.name, bak.name)

            write_ascii_tiff_comment(f, new_xml)
            logging.info("Updated XML for %s", f.name)
            successes += 1
        except Exception as exc:
            logging.exception("Failed to write TIFF comment for %s: %s", f, exc)

    logging.info(f"Completed {output_dir.name}: {successes}/{len(files)} files updated")
    return True


def main(argv=None):
    """_summary_
    Run batch rename of morphology_focus OME-TIFFs across all output directories.
    """  # noqa: D205
    parser = argparse.ArgumentParser(
        description=(
            "Batch rename morphology_focus OME-TIFFs across Xenium output directories"
        )
    )
    parser.add_argument(
        "-d",
        "--data-dir",
        required=True,
        help="Path to the directory containing output-* folders",
    )
    parser.add_argument(
        "-m",
        "--mapping",
        default="./xenium_mapping.md",
        help="Path to the mapping file (default: ./xenium_mapping.md)",
    )
    parser.add_argument(
        "--rename",
        action="store_true",
        help="Actually perform the renames (default is dry-run)",
    )
    parser.add_argument(
        "--backup", action="store_true", help="Create backup copies before renaming"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without making changes",
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Overwrite existing files"
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")

    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s: %(message)s",
    )

    data_dir = Path(args.data_dir)
    if not data_dir.exists():
        logging.error("Data directory does not exist: %s", data_dir)
        sys.exit(2)

    mapping_file = Path(args.mapping)
    if not mapping_file.exists():
        logging.error("Mapping file not found: %s", mapping_file)
        sys.exit(2)

    mapping = parse_mapping_file(mapping_file)
    logging.info(f"Loaded mapping with {len(mapping)} entries")

    # Find all output-* directories
    output_dirs = []
    for item in data_dir.iterdir():
        if item.is_dir() and item.name.startswith("output-"):
            output_dirs.append(item)

    if not output_dirs:
        logging.error("No output-* directories found in %s", data_dir)
        sys.exit(2)

    logging.info(f"Found {len(output_dirs)} output directories to process")

    # Process each output directory
    successful_dirs = 0
    for output_dir in sorted(output_dirs):
        if process_output_directory(output_dir, mapping, args):
            successful_dirs += 1

    logging.info(
        f"All {successful_dirs}/{len(output_dirs)} directories processed successfully"
    )


if __name__ == "__main__":
    main()
