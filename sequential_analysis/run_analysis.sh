
line_counter=1

# Loop through each line 
while IFS= read -r new_input_path; do
  echo "Processing line $line_counter with path: $new_input_path"

  #changes input in the betaconfig
  python3 change_inputfile_betaconfig.py "$line_counter" "run_list_100.txt"
  cd ../
  root -l -b -q Compilatore.C
  cd sequential_analysis/
  line_counter=$((line_counter + 1))

  echo "" 
done < "run_list_100.txt"

echo "Finished processing all lines in input_files.txt."
