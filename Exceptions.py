
class FileNotFoundError(Exception):
    def __init__(self):
        pass


class FileParseError(Exception):
    def __init__(self):
        pass


class SMPSCountError(Exception):
    def __init__(self):
        pass


class CCNCCountError(Exception):
    def __init__(self):
        pass


class NormalizedConcentrationError(Exception):
    def __init__(self):
        pass


class SigmoidFitGetParameterError(Exception):
    def __init__(self):
        pass


class SigmoidLineFitError(Exception):
    def __init__(self):
        pass


class ProgressBarInterruptException(Exception):
    def __init__(self):
        pass