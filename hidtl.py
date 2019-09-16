import sys
from src.hidtl_model import *


def default(str):
    return str + ' [Default: %default]'

def parseDistributionArgs(str):
    if str == None: return {}
    pieces = str.split(',')
    opts = {}
    for p in pieces:
        if '=' in p:
            key, val = p.split('=')
        opts[key] = int(val)
    return opts

def readCommand(argv):
    """
    Processes the command used to run HIDTLModel from the command line.
    """
    from optparse import OptionParser
    usageStr = """
    USAGE:      python hidtl.py <options>
    EXAMPLES:   (1) python hidtl.py
                    - runs a model
                (2) python hidtl.py --input data/species_tree.txt
                OR  python hidtl.py -i data/species_tree.txt
    """
    parser = OptionParser(usageStr, add_help_option=False)

    parser.add_option('--help', action='store_true', help='show this help message')   

    parser.add_option('-i', '--inputFile', dest='inputFile',
                      help='the path to an input file of a species tree', metavar='INPUT_FILE')
    parser.add_option('-c', '--coalescentArgs', dest='coalescentArgs',
                      help=default('the parameters of the gamma distribution for coalescent, '
                      'e.g., "shape=val1,scale=val2", or "const=val" if constant applies'),
                      default='shape=1,scale=1')
    parser.add_option('-d', '--duplicationArgs', dest='duplicationArgs',
                      help=default('the parameters of the gamma distribution for duplication event, '
                      'e.g., "shape=val1,scale=val2", or "const=val" if constant applies'),
                      default='shape=1,scale=1')
    parser.add_option('-t', '--transferArgs', dest='transferArgs',
                      help=default('the parameters of the gamma distribution for transfer event, '
                      'e.g., "shape=val1,scale=val2", or "const=val" if constant applies'),
                      default='shape=1,scale=1')
    parser.add_option('-l', '--lossArgs', dest='lossArgs',
                      help=default('the parameters of the gamma distribution for loss event, '
                      'e.g., "shape=val1,scale=val2", or "const=val" if constant applies'),
                      default='shape=1,scale=1') 
    parser.add_option('-h', '--hemiplasy', type='int', dest='hemiplasy',
                      help=default('hemiplasy option, 0 or 1'), metavar='HEMIPLASY',
                      default=1)      
    parser.add_option('-r', '--recombination', type='int', dest='recombination',
                      help=default('recombination option, 0 or 1'), metavar='RECOMBINATION',
                      default=1)   

    options, otherjunk = parser.parse_args(argv)
    if len(otherjunk) != 0:
        raise Exception('Command line input not understood: ' + str(otherjunk))
    args = dict()

    if options.help:
        parser.print_help()
        sys.exit()

    # input file (a species tree in newick format)
    if not options.inputFile:
        parser.error('The input filename not given')
    args['inputFile'] = options.inputFile

    # distribution arguments
    args['coalescentArgs'] = parseDistributionArgs(options.coalescentArgs)
    args['duplicationArgs'] = parseDistributionArgs(options.duplicationArgs)
    args['transferArgs'] = parseDistributionArgs(options.transferArgs)
    args['lossArgs'] = parseDistributionArgs(options.lossArgs)

    # hemiplasy option
    if options.hemiplasy != 0 and options.hemiplasy != 1:
        parser.error('Invalid hemiplasy option: ' + str(options.hemiplasy))
    args['hemiplasy'] = True if options.hemiplasy == 1 else False

    # recombination option
    if options.recombination != 0 and options.recombination != 1:
        parser.error('Invalid recombination option: ' + str(options.recombination))
    args['recombination'] = True if options.recombination == 1 else False

    return args

def runModel(inputFile, coalescentArgs, duplicationArgs, transferArgs, lossArgs,
             hemiplasy, recombination):
    model = HIDTLModel()

    model.setParameters(coalescent=coalescentArgs, 
                        duplication=duplicationArgs, 
                        transfer=transferArgs, 
                        loss=lossArgs,
                        hemiplasy=hemiplasy,
                        recombination=recombination)

    model.readSpeciesTree(inputFile)

    model.createHaplotypeTree()

if __name__ == '__main__':
    """
    The hidtl function called when hidtl.py is run
    from the command line:

    > python hidtl.py

    See the usage string for more details.

    > python hidtl.py --help
    """
    args = readCommand(sys.argv[1:]) # Get game components based on input
    runModel(**args)

    pass
