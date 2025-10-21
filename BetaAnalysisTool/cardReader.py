# cardReader.py

import re

def read_text_card(file_path):
  config = {}

  channel_type_mapping = {
      "DUT": 1,
      "MCP": 2,
      "REF": 3
  }
  channel_area_to_charge_mapping = {
        "SC": 4.7,
        "Mig": 5
  }
  channels = [[0, 1]] * 8

  plot_flags = {
    "tmax": False,
    "pmax": False,
    "negpmax": False,
    "amplitude": False,
    "risetime": False,
    "charge": False,
    "rms": False,
    "timeres": False,
    "discretisation": False,
    "waveform": False,
    "all_with_no_plots": False,
  }

  plot_params = {
    "tmax_params": None,
    "pmax_params": None,
    "negpmax_params": None,
    "risetime_params": None,
    "charge_params": None,
    "rms_params": None,
    "timeres_params": None
  }

  MCP_specs = None

  with open(file_path, 'r') as f:
    current_key = None  # Track the current key being processed
    current_value = []  # Collect multi-line values
    thickness_info = []
    for line in f:
      line = line.strip()
      if not line or line.startswith('#'):  # Skip empty lines and comments
        continue

      match = re.match(r'^(\w+)\s*=\s*(.+)$', line)
      if match:
        if current_key and current_value:
          if current_key == "files":
            config[current_key] = "".join(current_value).strip('",').split(',')
          else:
            config[current_key] = "".join(current_value).strip()
          current_key = None
          current_value = []

        key, value = match.groups()
        value = value.strip().strip('"').strip("'")

        if key == "files":  # Handle multi-line `files`
          current_key = key
          current_value.append(value)

        elif key == "run_safe_mode":
          safemode = value.lower() == "true"

        elif key.startswith("CH_") and key[3:].isdigit():  # Handle CH_ keys
          index = int(key[3:]) - 1  # Convert to 0-based index
          parts = [part.strip() for part in value.split(',')]
          type_str = parts[0]
          additional_str = parts[1] if len(parts) > 1 else ""
          thickness_str = parts[2] if (len(parts) > 2) & (type_str.upper() != "MCP") else "nDUT"

          channel_type = channel_type_mapping.get(type_str.upper(), 0)
          channel_value = channel_area_to_charge_mapping.get(additional_str, 1)
          thickness_info.append(thickness_str)

          if (type_str.upper() == "MCP"):
            if (parts[2] == 0) & (parts[3] == 0):
              MCP_specs = (0, 0)
            else:
              MCP_specs = (float(parts[2]), float(parts[3]))

          channels[index] = [channel_type, channel_value, None]

        elif key.startswith("CH") and key.endswith("_cut"):
          channel_index = int(key[2]) - 1

          match = re.match(r"^\s*(\[\s*(?:-?\d+(?:\.\d+)?\s*,\s*)*-?\d+(?:\.\d+)?\s*\]|\[\s*\]|0)\s*,\s*(.*)$", value)
          if not match:
            raise ValueError(f"Invalid format for {key}: Must start with an array '[x,y,...]', '[]', or '0'.")

          raw_lower_bound = match.group(1).strip()
          remaining_values = match.group(2).strip()

          if raw_lower_bound == "[]" or raw_lower_bound == "0":
            lower_bound = []
          else:
            lower_bound = list(map(float, raw_lower_bound.strip("[]").split(",")))

            if len(lower_bound) != len(config.get('files', [])):
              raise ValueError(
                  f"Invalid length for lower_bound in {key}: Must match the number of files ({len(config['files'])})."
              )

          remaining_parts = remaining_values.split(",")
          if len(remaining_parts) != 4:
            raise ValueError(f"Invalid format for {key}: Must contain exactly 4 additional comma-separated values.")

          try:
            upper_bound, additional_condition, tlow, thigh = map(float, remaining_parts)
          except ValueError:
            raise ValueError(f"Invalid format for {key}: The last 4 values must all be floats.")

          channels[channel_index][2] = (lower_bound, upper_bound, additional_condition, tlow, thigh)

        elif key in plot_flags:  # Handle plot flags
          plot_flags[key] = value.lower() == "true"
        elif key.endswith("_nB_xL_xU"):  # Handle plot parameters
          param_key = key.split("_nB_xL_xU")[0]
          if plot_flags.get(param_key, False) == True:  # Check if this plot is enabled
            nBins, xLower, xUpper = map(float, value.split(","))
            plot_params[param_key+"_params"] = (int(nBins), xLower, xUpper)
          elif (plot_flags.get("amplitude", False)) and (param_key == "pmax"):
            nBins, xLower, xUpper = map(float, value.split(","))
            plot_params["pmax_params"] = (int(nBins), xLower, xUpper)
        else:  # Handle generic key-value pairs
          config[key] = value
      elif current_key:  # Handle continuation lines
        current_value.append(line.strip())

    if current_key and current_value:
        if current_key == "files":
          config[current_key] = "".join(current_value).strip('",').split(',')
        else:
          config[current_key] = "".join(current_value).strip()

  config['channels'] = channels
  config.update(plot_flags)
  config.update(plot_params)
  return config, thickness_info, MCP_specs, safemode
