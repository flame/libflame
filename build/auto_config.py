"""Copyright (C) 2023, Advanced Micro Devices, Inc. All Rights Reserved"""

import platform
import re
import subprocess
import sys

def config_check():
    # Execute wmic shell command with sub-process
    global model, family, vendor, stepping
    if 'Windows' in platform.system():
        result = subprocess.Popen(
            'wmic cpu get caption', shell=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        result = result[0].decode('utf-8')
    # Replace the newline character with empty char
        result = result.replace('\n', '')

    # parse the string into list of string
        parse_string = result.split(" ")

    # Strip the empty strings from list
        parse_string = [data for data in parse_string if data.strip()]

        vendor = parse_string[1]
        family = int(parse_string[3])
        model = int(parse_string[5])
        stepping = int(parse_string[7])

    elif 'Linux' in platform.system():
        result = subprocess.Popen(
            'lscpu', shell=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        stepping = int(re.findall(r'\WStepping:.*', result[0].decode(
            'utf-8'), re.MULTILINE)[0].strip('\n').split(' ')[-1])
        family = int(re.findall(r'\WCPU family:.*', result[0].decode(
            'utf-8'), re.MULTILINE)[0].strip('\n').split(' ')[-1])
        model = int(re.findall(r'\WModel:.*', result[0].decode('utf-8'),
                               re.MULTILINE)[0].strip('\n').split(' ')[-1])
        vendor = re.findall(r'\WModel name:.*', result[0].decode(
            'utf-8'), re.MULTILINE)[0].strip('\n').split('Model name:')[-1]
    # AMD family numbers
    # Zen/Zen+/Zen2 (0x17) and Zen3/Zen4 (0x19) family numbers
    zen_family = [23, 25]
    # Bulldozer / Piledriver / Steamroller / Excavator family number
    amd_family = 21

    # AMD CPUID model numbers
    zen_model = [48, 255]
    zen2_model = [0, 255]
    zen3_model = [(0, 15), (32, 95)]
    zen4_model = [(16, 31), (96, 175)]
    excavator_model = [96, 127]
    steamroller_model = [48, 63]
    piledriver_model = [2, 16, 31]
    bulldozer_model = [0, 1]

    # Check the CPU configuration Intel64/AMD64
    if vendor.count("Intel64"):
        return
    elif 'AMD' in vendor:  # .count("AMD64"):
        # Check the AMD family name
        if family == zen_family[0]:
            if zen_model[0] <= model <= zen_model[1]:
                family="zen2"
            elif zen2_model[0] <= model <= zen2_model[1]:
                family="zen"
            else:
                print("Unknown model number")
        elif family == zen_family[1]:
            if (zen3_model[0][0] <= model <= zen3_model[0][1]) or (
                    zen3_model[1][0] <= model <= zen3_model[1][1]):
                family="zen3"
            elif (zen4_model[0][0] <= model <= zen4_model[0][1]) or (
                    zen4_model[1][0] <= model <= zen4_model[1][1]):
                family="zen4"
            else:
                print("Unknown model number zen4")
        elif family == amd_family:
            # check for specific models of excavator family
            if excavator_model[0] <= model <= excavator_model[1]:
                family="excavator"
            # check for specific models of steamroller family
            elif steamroller_model[0] <= model <= steamroller_model[1]:
                family="steamroller"
            # check for specific models of piledriver family
            elif model == piledriver_model[0] or (
                    piledriver_model[1] <= model <= piledriver_model[2]):
                family="piledriver"
            # check for specific models of bulldozer family
            elif model == bulldozer_model[0] or model == bulldozer_model[1]:
                family="bulldozer"
            else:
                print("Unknown model number")
        else:
            print("Unknown family")
    else:
        print("UNKNOWN CPU")
    return family

# Function call for config family names
config = config_check()
print(config)
