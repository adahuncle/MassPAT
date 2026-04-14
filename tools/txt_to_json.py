import json
import re


def extract_numeric(value):
    """
    take '37.97631922(21)' -> 37.97631922
    """
    if not value:
        return None
    match = re.match(r"([-+]?\d*\.\d+|\d+)", value)
    return float(match.group(0)) if match else None


def parse_value_with_uncertainty(value):
    """
    Parse values like '30.975363194(46)' into:
      value: 30.975363194
      uncertainty: 0.000000046

    Parenthetical digits are interpreted as uncertainty in the last digits.
    """
    if not value:
        return {
            "value": None,
            "uncertainty": None,
            "value_text": None,
            "uncertainty_text": None,
            "plus_minus": None,
        }

    text = value.strip()
    match = re.match(r"^([-+]?\d+(?:\.\d+)?)(?:\((\d+)\))?$", text)
    if not match:
        numeric = extract_numeric(text)
        return {
            "value": numeric,
            "uncertainty": None,
            "value_text": str(numeric) if numeric is not None else None,
            "uncertainty_text": None,
            "plus_minus": None,
        }

    base_text = match.group(1)
    unc_digits = match.group(2)
    base_value = float(base_text)

    if unc_digits is None:
        return {
            "value": base_value,
            "uncertainty": None,
            "value_text": base_text,
            "uncertainty_text": None,
            "plus_minus": None,
        }

    decimals = len(base_text.split(".")[1]) if "." in base_text else 0
    uncertainty = int(unc_digits) * (10 ** (-decimals))
    if decimals > 0:
        uncertainty_text = f"{uncertainty:.{decimals}f}"
    else:
        uncertainty_text = str(int(uncertainty))

    return {
        "value": base_value,
        "uncertainty": uncertainty,
        "value_text": base_text,
        "uncertainty_text": uncertainty_text,
        "plus_minus": f"{base_text}±{uncertainty_text}",
    }


def parse_standard_atomic_weight(value):
    """
    Parse standard atomic weight into either:
      - range: [min,max]
      - specific: value(unc)
    """
    if not value:
        return {
            "kind": None,
            "range_min": None,
            "range_max": None,
            "specific_value": None,
            "uncertainty": None,
            "plus_minus": None,
        }

    text = value.strip()
    range_match = re.match(r"^\[\s*([-+]?\d+(?:\.\d+)?)\s*,\s*([-+]?\d+(?:\.\d+)?)\s*\]$", text)
    if range_match:
        return {
            "kind": "range",
            "range_min": float(range_match.group(1)),
            "range_max": float(range_match.group(2)),
            "specific_value": None,
            "uncertainty": None,
            "plus_minus": None,
        }

    parsed_specific = parse_value_with_uncertainty(text)
    return {
        "kind": "specific" if parsed_specific["value"] is not None else None,
        "range_min": None,
        "range_max": None,
        "specific_value": parsed_specific["value"],
        "uncertainty": parsed_specific["uncertainty"],
        "plus_minus": parsed_specific["plus_minus"],
    }


def parse_nist_ascii(file_path):
    with open(file_path, "r", encoding="utf-8") as f:
        text = f.read()

    # split blocks by blank lines
    blocks = re.split(r"\n\s*\n", text.strip())

    records = []

    for block in blocks:
        row = {
            "atomic_number": None,
            "atomic_symbol": "",
            "mass_number": None,
            "relative_atomic_mass": "",
            "relative_atomic_mass_value": None,
            "relative_atomic_mass_uncertainty": None,
            "relative_atomic_mass_plus_minus": None,
            "isotopic_composition": "",
            "isotopic_composition_value": None,
            "standard_atomic_weight": "",
            "standard_atomic_weight_value": None,
            "standard_atomic_weight_kind": None,
            "standard_atomic_weight_range_min": None,
            "standard_atomic_weight_range_max": None,
            "standard_atomic_weight_specific_value": None,
            "standard_atomic_weight_uncertainty": None,
            "standard_atomic_weight_plus_minus": None,
            "notes": ""
        }

        for line in block.splitlines():
            if "=" not in line:
                continue

            key, value = line.split("=", 1)
            key = key.strip()
            value = value.strip()

            if key == "Atomic Number":
                row["atomic_number"] = int(value) if value else None

            elif key == "Atomic Symbol":
                row["atomic_symbol"] = value

            elif key == "Mass Number":
                row["mass_number"] = int(value) if value else None

            elif key == "Relative Atomic Mass":
                row["relative_atomic_mass"] = value
                parsed_mass = parse_value_with_uncertainty(value)
                row["relative_atomic_mass_value"] = parsed_mass["value"]
                row["relative_atomic_mass_uncertainty"] = parsed_mass["uncertainty"]
                row["relative_atomic_mass_plus_minus"] = parsed_mass["plus_minus"]

            elif key == "Isotopic Composition":
                row["isotopic_composition"] = value
                row["isotopic_composition_value"] = extract_numeric(value)

            elif key == "Standard Atomic Weight":
                row["standard_atomic_weight"] = value
                row["standard_atomic_weight_value"] = extract_numeric(value)
                parsed_weight = parse_standard_atomic_weight(value)
                row["standard_atomic_weight_kind"] = parsed_weight["kind"]
                row["standard_atomic_weight_range_min"] = parsed_weight["range_min"]
                row["standard_atomic_weight_range_max"] = parsed_weight["range_max"]
                row["standard_atomic_weight_specific_value"] = parsed_weight["specific_value"]
                row["standard_atomic_weight_uncertainty"] = parsed_weight["uncertainty"]
                row["standard_atomic_weight_plus_minus"] = parsed_weight["plus_minus"]

            elif key == "Notes":
                row["notes"] = value

        # only keep valid isotope rows
        if row["atomic_symbol"] and row["mass_number"] is not None:
            row["isotope_label"] = f'{row["atomic_symbol"]}-{row["mass_number"]}'
            records.append(row)

    return records


def save_json(records, output_path="isotopes.json"):
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(records, f, indent=2)


if __name__ == "__main__":
    from pathlib import Path
    _root = Path(__file__).resolve().parent.parent
    input_file = str(_root / "data" / "reference" / "isotopic_nist_info.txt")
    output_file = str(_root / "data" / "reference" / "isotopes.json")

    data = parse_nist_ascii(input_file)
    save_json(data, output_file)

    print(f"done. saved {len(data)} isotopes to {output_file}")