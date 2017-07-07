
class FileNotFoundError(Exception):
    def __init__(self):
        pass


class FileProcessingError(Exception):
    def __init__(self):
        pass


class DNlogDataError(Exception):
    def __init__(self):
        pass


class DataError(Exception):
    def __init__(self):
        pass


class DataMatchingError(Exception):
    def __init__(self):
        pass


class DataPreparationError(Exception):
    def __init__(self):
        pass


class ScanDataError(Exception):
    def __init__(self):
        pass


class SigmoidFitError(Exception):
    def __init__(self):
        pass


class ProgressBarInterruptException(Exception):
    def __init__(self):
        pass