def delete_last_n_lines(filename, n):
    """
    Delete the last n lines from a given file.

    Args:
        filename (str): The name of the file.
        n (int): The number of lines to delete.

    Returns:
        None
    """
    try:
        # Read all lines from the file
        with open(filename, 'r') as file:
            lines = file.readlines()

        # Ensure n is not greater than the number of lines in the file
        n = min(n, len(lines))

        # Write back all lines excluding the last n lines
        with open(filename, 'w') as file:
            file.writelines(lines[:-n])

        print(f"Last {n} lines deleted from {filename}.")
    except FileNotFoundError:
        print(f"File '{filename}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

# Example usage:
# Delete the last 2 lines from a file named 'example.txt'
delete_last_n_lines('data/expr_sims/20240114_readphase_free_hemisphere/ciliate_639fil_40961blob_15.00R_0.042torsion_blob_forces.dat', 95)
