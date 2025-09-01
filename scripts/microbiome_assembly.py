#!/usr/bin/env python3
import sys
import os
import re
import subprocess

def main():
    """
    Main function to process the TSV file and execute generated commands.
    """
    # === Argument Handling ===
    # The script requires exactly two command-line arguments.
    if len(sys.argv) != 3:
        # Print usage information to standard error and exit.
        print(f"Usage: {sys.argv[0]} <input_tsv_file> <output_directory>", file=sys.stderr)
        sys.exit(1)

    input_file = sys.argv[1]
    output_root = sys.argv[2]

    # Check if the input file exists.
    if not os.path.isfile(input_file):
        print(f"Error: Input file not found: {input_file}", file=sys.stderr)
        sys.exit(1)

    # === Directory Management ===
    # Create the output directory if it doesn't exist.
    os.makedirs(output_root, exist_ok=True)

    # Change the current working directory to the output directory.
    # All subsequent file operations will be relative to this path.
    try:
        os.chdir(output_root)
        print(f"? Output Directory: {os.getcwd()}")
    except OSError as e:
        print(f"Error: chdir to {output_root} failed: {e}", file=sys.stderr)
        sys.exit(1)

    # === Command Generation ===
    # This list will store all the commands to be executed.
    commands_to_run = []
    
    try:
        # Open the input TSV file for reading.
        with open(input_file, 'r', encoding='utf-8') as f_in:
            # Read the header line.
            header_line = f_in.readline()
            if not header_line:
                print("Warning: Input file is empty or header is missing.", file=sys.stderr)
                return
            
            headers = header_line.strip().split('\t')

            # Process each data line in the TSV.
            for line in f_in:
                cols = line.strip().split('\t')
                if len(cols) < len(headers):
                    continue # Skip malformed lines

                # Assume the last column is the taxid.
                taxid = cols[-1]
                
                # Skip if taxid is not a valid number or is 'Unassigned'.
                if not taxid or not taxid.isdigit() or taxid == 'Unassigned':
                    continue

                # Find the lowest taxonomic rank name that is not 'Unassigned'.
                name, header = '', ''
                # Iterate backwards from column index 8 to 1.
                for i in range(min(8, len(cols) - 1), 0, -1):
                    # Check if the column exists and has a valid value.
                    if cols[i] and cols[i].strip() and cols[i] != 'Unassigned':
                        name = cols[i]
                        header = headers[i]
                        break # Found a valid name, exit the loop.

                # If no valid name was found, skip this line.
                if not name:
                    continue
                
                # Sanitize the name and header to be used in file/directory names.
                # Replace forbidden characters and whitespace with underscores.
                safe_name = re.sub(r'[\\/:*?"<>|]', '_', name)
                safe_name = re.sub(r'\s+', '_', safe_name)
                safe_header = re.sub(r'[\\/:*?"<>|]', '_', header)
                safe_header = re.sub(r'\s+', '_', safe_header)

                # Create subdirectory and define the final output file path.
                dirname = f"{safe_header}_{safe_name}"
                os.makedirs(dirname, exist_ok=True)
                output_path = os.path.join(dirname, f"{safe_header}_{safe_name}.tsv")

                # Generate the shell command.
                cmd = f"datasets summary genome taxon {taxid} --as-json-lines --reference | dataformat tsv genome --fields accession,organism-name,assminfo-level,organism-tax-id > {output_path}"
                commands_to_run.append(cmd)

    except FileNotFoundError:
        print(f"Error: Could not open input file {input_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred during file processing: {e}", file=sys.stderr)
        sys.exit(1)

    # Write all generated commands to 'commands.txt'.
    try:
        with open('commands.txt', 'w', encoding='utf-8') as f_cmd:
            for cmd in commands_to_run:
                f_cmd.write(cmd + '\n')
    except IOError as e:
        print(f"Error: Cannot create commands.txt: {e}", file=sys.stderr)
        sys.exit(1)


    # === Command Execution ===
    print("--- Executing generated commands ---")
    for cmd in commands_to_run:
        if not cmd.strip():
            continue
        
        print(f"Do: {cmd}")
        # Execute the command using the shell. `shell=True` is required for redirection (`>`).
        # This mimics the behavior of Perl's `system()`.
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        # If the command fails, print a warning to stderr but continue execution.
        if result.returncode != 0:
            print(f"Warning: command failed (exit code {result.returncode}): {cmd}", file=sys.stderr)
            if result.stderr:
                print(f"Stderr: {result.stderr.strip()}", file=sys.stderr)

    print("--- Done ---")


if __name__ == "__main__":
    main()