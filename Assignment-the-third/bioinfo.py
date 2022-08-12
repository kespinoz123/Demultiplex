# Author: Itzel Espinoza kespinoz@uoregon.edu


def convert_phred(letter: str) -> int:
    """Converts a single character into a phred score"""
    convert_phred=(ord(letter)-33) 
    return convert_phred

