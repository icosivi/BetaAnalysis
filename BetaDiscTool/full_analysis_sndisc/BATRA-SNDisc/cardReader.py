# cardReader.py

import re
import ast

def read_text_card(file_path):
  config = {}

  channel_type_mapping = {
      "DUT": 1,
      "MCP": 2,
      "REF": 3
  }
  channel_area_to_charge_mapping = {
        "SC": 4.7,
        "Mig": 5,
        "None": 1
  }
  channels = [[0, 1]] * 8

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

        elif key.startswith("CH_"):  # Handle CH_ keys
          index = int(key[3:]) - 1  # Convert to 0-based index
          pattern = r',\s*(?![^\[\]()]*[\]\)])'
          parts = re.split(pattern, value)
          channel_type = channel_type_mapping.get(parts[0].upper(), 0)
          
          if (channel_type == 1):
            if len(parts) != 4:
              raise ValueError(
                f"Invalid length for setting in {key}: Must specify sensor type, board, thickness, and array of ansatz PMAX values."
              )
            channel_value = channel_area_to_charge_mapping.get(parts[1], 1)
            thickness_str = parts[2]
            channel_selections = ast.literal_eval(parts[3])
            if len(channel_selections) != len(config.get('files', [])):
              raise ValueError(
                f"Invalid length for lower_bound in {key}: Must match the number of files ({len(config['files'])})."
              )
          elif (channel_type == 2):
            if len(parts) != 4:
              raise ValueError(
                f"Invalid length for setting in {key}: Must specify sensor type, board, MCP time resolution and uncertainty (as a pair), and min/max PMAX range (as a pair)."
              )
            channel_value = channel_area_to_charge_mapping.get(parts[1], 1)
            MCP_specs = ast.literal_eval(parts[2])
            channel_selections = ast.literal_eval(parts[3])
            thickness_str = "nDUT"
          elif (channel_type == 3):
            channel_value = channel_area_to_charge_mapping.get(parts[1], 1)
            thickness_str = "nDUT"
            channel_selections = []
          else:
            additional_str = ""
            thickness_str = "nDUT"
            channel_selections = []
            
          thickness_info.append(thickness_str)

          channels[index] = [channel_type, channel_value, channel_selections]
          
        elif key.startswith("first_pass"):
          pass_flag_and_pmax = [part.strip() for part in value.split(',')]
          pass_flag_and_pmax[0] = pass_flag_and_pmax[0].lower() == "true"
          pass_flag_and_pmax[1] = int(pass_flag_and_pmax[1])
        
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
  config['pass_criteria'] = pass_flag_and_pmax
  return config, thickness_info, MCP_specs
