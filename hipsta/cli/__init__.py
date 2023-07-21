"""
hipsta CLI
"""

from ..hipsta import _parse_arguments, run


def main():

    # message

    print('')
    print('----------------------------------------')
    print('Hippocampal shape and thickness analysis')
    print('----------------------------------------')
    print('')

    # parse arguments

    args = _parse_arguments()

    # return help or run analysis

    if args.more_help is True:
        get_help()
    else:
        run(args)