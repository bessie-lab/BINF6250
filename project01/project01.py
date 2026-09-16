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

def read_file(file):
    
    """
    Reads a vcf file, takes each line, runs it through the parse line function above and updates a growing dictionary 
    with the occurrence of a disease as it shows up in the vcf. 
    
    Returns a dictionary with all the times a disease is detected in the vcf file based on parameters from parse_line()

    If the AF_EXAC parameter isn't found or is greater than or equal to 0.0001, the disease isn't counted. 
    If the AF_EXAC parameter is below 0.0001, the variant is considered rare and all associated diseases are reported.
    """
    
    vcf_disease_count = {}

    with open(file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            else:
                parse_outputs = parse_line(line)

                # If a disease was reported from the parse_line function in a list
                if parse_outputs is not None:
                    for disease in parse_outputs:
                        if disease not in vcf_disease_count.keys():
                            vcf_disease_count[disease] = 1
                        else:
                            vcf_disease_count[disease] += 1

    return vcf_disease_count

if __name__ == "__main__":
    pprint(read_file("clinvar_20190923_short.vcf"))
