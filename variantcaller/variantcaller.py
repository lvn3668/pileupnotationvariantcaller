def checkforreferenceallelecall(uniquebasescalled: list[str]) -> bool:
    """

    :param uniquebasescalled:
    :return:
    """
    try:
        if len(uniquebasescalled) == 1 and uniquebasescalled[0] == 0.:
            return True
        else:
            return False
    except Exception as exception:
        print(type(exception))  # the exception instance
        print(exception.args)  # arguments stored in .args
        print(exception)  # __str__ allows args to be printed directly,
        x, y = exception.args  # unpack args
        print('x =', x)
        print('y =', y)


def checkforaltallelecallinpileup(uniquebasescalled: list[str]) -> bool:
    """

    :param uniquebasescalled:
    :return:
    """
    try:
        if len(uniquebasescalled) == 1 and uniquebasescalled[0] == 1.:
            return True
        else:
            return False
    except Exception as exception:
        print(type(exception))  # the exception instance
        print(exception.args)  # arguments stored in .args
        print(exception)  # __str__ allows args to be printed directly,
        x, y = exception.args  # unpack args
        print('x =', x)
        print('y =', y)


def checkforalthomozygouscall(uniquebasescall: list[str], frequencies) -> bool:
    """

    :param uniquebasescall:
    :param frequencies:
    :return:
    """
    try:
        if (len(uniquebasescall) == 2 and
                ((frequencies[1][1] / (frequencies[0][1] + frequencies[1][1])) > 0.75)):
            return True
        else:
            return False
    except Exception as exception:
        print(type(exception))  # the exception instance
        print(exception.args)  # arguments stored in .args
        print(exception)  # __str__ allows args to be printed directly,
        x, y = exception.args  # unpack args
        print('x =', x)
        print('y =', y)


def checkifATGClow(uniquebasescalled: list[str], frequencies: object) -> bool:
    """
    :return:
    :param uniquebasescalled:
    :param frequencies:
    :return:
    """
    try:
        if (len(uniquebasescalled) == 2 and
                ((frequencies[1][1] / (frequencies[0][1] + frequencies[1][1])) < 0.25)):
            return True
        else:
            return False
    except Exception as exception:
        print(type(exception))  # the exception instance
        print(exception.args)  # arguments stored in .args
        print(exception)  # __str__ allows args to be printed directly,
        x, y = exception.args  # unpack args
        print('x =', x)
        print('y =', y)


def checkifATGCHeterozygous(uniquebasescalled: list[str], frequencies) -> bool:
    """

    :rtype: bool
    :param uniquebasescalled:
    :param frequencies:
    :return:
    """
    try:
        if (len(uniquebasescalled) == 2 and
                ((frequencies[1][1] / (frequencies[0][1] + frequencies[1][1])) >= 0.25)
                &
                ((frequencies[1][1] / (frequencies[0][1] + frequencies[1][1])) <= 0.75)):
            return True
        else:
            return False
    except Exception as exception:
        print(type(exception))  # the exception instance
        print(exception.args)  # arguments stored in .args
        print(exception)  # __str__ allows args to be printed directly,
        x, y = exception.args  # unpack args
        print('x =', x)
        print('y =', y)
