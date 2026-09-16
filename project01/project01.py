#!/usr/bin/env python
from pprint import pprint
EXCLUDED_DISEASES = {"not_specified", "not_provided"}
INFO_COLUMN = 7


def parse_line(line):
    """
    Parse a single VCF data line .It takes a string as an argument and returns
        list: A list of disease names (str) associated with the variant
              if the variant is rare (AF_EXAC < 0.0001).
              An empty list if the variant is NOT rare.
              None if AF_EXAC is not present in the INFO field (skip).
    """
    fields = line.rstrip("\n").split("\t")

    if len(fields) <= INFO_COLUMN:
        return None

    info_field = fields[INFO_COLUMN]

    info_dict = {}
    for entry in info_field.split(";"):
        if "=" in entry:
            key, value = entry.split("=", 1)
            info_dict[key] = value


    if "AF_EXAC" not in info_dict:
        return None

    af_exac = float(info_dict["AF_EXAC"])

    if af_exac < 0.0001:
        clndn_raw = info_dict.get("CLNDN", "")
        diseases = clndn_raw.split("|") if clndn_raw else []
        diseases = [d for d in diseases if d not in EXCLUDED_DISEASES]
        return diseases


    return []
    
