"""
Utility Classes
"""


class Value:
    """ Container for value and unit  """

    def __init__(self, value: float, unit='') -> None:
        self.value = value
        self.unit = unit


class FileUrl:
    """ Container for local file location and URL reference  """

    def __init__(self, file, url=None) -> None:
        self.file = file
        self.url = url
