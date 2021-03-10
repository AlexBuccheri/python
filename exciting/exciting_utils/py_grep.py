"""
Python wrapper for grep and other useful spec
"""
import subprocess


# TODO Fix to take **args too and refactor line_number
def grep(string: str, fname: str, **kwargs):
    """
    Grep wrapper

    :param string:
    :param fname:
    :param kwargs:
    :return:
    """
    grep_options = {'n_lines_before': '-B',
                    'n_lines_after': '-A',
                    'line_number': '-n'}

    opts = ''
    for key, value in kwargs.items():
        opts += grep_options[key] + ' ' + str(value) + ' '

    grep_str = "grep " + opts + " '" + string + "' " + fname

    try:
        output = subprocess.check_output(grep_str, shell=True).decode("utf-8")
    except subprocess.CalledProcessError as grepexc:
        print("subprocess error:", grepexc.returncode, "grep found:", grepexc.output)
        #TODO throw error if string not returned
        output = grepexc.output

    return output


#
#
# def number_of_string_matches(string, fname):
#     grep_str = "grep -o '" + string + "' "+fname+" | wc -l"
#     n_matches = subprocess.check_output(grep_str, shell=True).decode("utf-8")
#     return int(n_matches)