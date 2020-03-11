import subprocess


# Python wrapper for grep
# Originally written in Work_2019/DL_POLY_DFTB+/python_plotting/extract.py
#
# Example usage:
# extract(">> Charges saved for restart in charges.bin", fname, n_lines_before=1)

grep_options = {'n_lines_before': '-B',
                'n_lines_after': '-A'}


def extract(string, fname, **kwargs):

    opts = ''
    for key, value in kwargs.items():
        opts += grep_options[key] + ' ' + str(value) + ' '

    grep_str = "grep " + opts + " '" + string + "' " + fname

    try:
        output = subprocess.check_output(grep_str, shell=True).decode("utf-8")
    except subprocess.CalledProcessError as grepexc:
        print("subprocess error:", grepexc.returncode, "grep found:", grepexc.output)
        output = grepexc.output

    return output


def number_of_string_matches(string, fname):
    grep_str = "grep -o '" + string + "' "+fname+" | wc -l"
    n_matches = subprocess.check_output(grep_str, shell=True).decode("utf-8")
    return int(n_matches)