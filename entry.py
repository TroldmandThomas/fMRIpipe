import sys, os
import argparse
import numpy as np
import pandas as pd
import pipeline.loadmatrix as lm 
import pipeline.obtain_estimates as oe
import statistics.get_ttest as gtt
import statistics.draw_graphs as dg
import statistics.glm as glm

#initialize parser
parser = argparse.ArgumentParser()

#initialize the positional argument "mode"
parser.add_argument("mode", help="Choose either full, estimate, ttest, plots or glm.")

#initialize the optional arguments
parser.add_argument('-mat', nargs='?', help="The MATLAB Conn file containing the matrices.")
parser.add_argument('-id', nargs='?', help="The CSV file containing the ID for the subjects.")
parser.add_argument('-thr', nargs='?', help="The threshold range or list of thresholds.")
parser.add_argument('-cut', nargs='?', 
         help="The part of the matrix that needs to extracted, default is the full matrix. ")
parser.add_argument('-ws', nargs='?', help="'W' for winter, 'S' for summer.")
parser.add_argument('-dir', nargs='?', help="Path to the estimate files.")
parser.add_argument('-out', nargs='?', help="Path to where the resulting CSV files should be written to. ")

args = parser.parse_args()


#blockPrint and enablePrint by courtesy of 
#https://stackoverflow.com/questions/8391411/suppress-calls-to-print-python
#FakeRainBrigand

# Disable printing function
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

# Restore printing function
def enablePrint():
    sys.stdout = sys.__stdout__

#to estimate the graph theory measures from the given matrices
def run_graph_estimates(pm,out=args.out):
    thresh_tok = args.thr.split(':')
    start = int(thresh_tok[0])
    end = int(thresh_tok[1])
    stride = int(thresh_tok[2])
    thresh_list = []
    for i in range(start,(end+stride),stride):
            thresh_list.append(float(i)/100)


    #print('Converting MATLAB matrices to NumPy arrays..')
    print('Running the graph theory estimations..')
    for thresh in thresh_list:
        print('Now processing threshold: ' +str(round(100 * thresh,2)) + '%')
        oe.obtain_estimates(pm, args.id, thresh, out)
    print("Graph theory estimates completed on all thresholds")


#pull out the MATLAB matrices from the Conn MATLAB file
def extract_matlab_mats():
    if not args.cut:
        size = 'full'
    else:
        size = args.cut
    cms = list(map(str, args.mat.strip('[]').split(',')))
    pm = lm.conn_interface(cms, size)

    return pm


#standard error message for when user forgets some parameter, or incorrectly entered
def error_msg():
    print(' ')
    print('**Something went wrong, please check input were correctly formatted,**')
    print('**or check that some parameter is not missing.**')
    print(' ')
    parser.print_help()



#running the full pipeline
if args.mode == 'full':

    try:
        pm = extract_matlab_mats()
        run_graph_estimates(pm, out=args.out)
        #get_ttest is called through draw_graphs
        direc = args.out + '/auto_results/'
        print('Drawing graphs..')
        dg.execute(path=direc, go=args.out, dest=args.out)

        print('Full pipeline run completed.')
    except:
        error_msg()


#running only the graph theory estimates on a MATLAB matrix
elif args.mode == 'estimate':

    try:
        pm = extract_matlab_mats()
        run_graph_estimates(pm)
        print('Done.')
    except:
        error_msg()

#running only the t-tests
elif args.mode == 'ttest':

    print('Performing t-tests..')
    if args.ws:
        gtt.gtt_main(WS=args.ws,path=args.dir, dest=args.out)
    else:
        gtt.gtt_main(path=args.dir, dest=args.out)
    print('Done.')

#running only the drawing of graphs (requires t-test to be run also)
elif args.mode == 'plots':

    print('Drawing plots..')
    dg.execute(path=args.dir, go=args.out)
    print('Done.')

elif args.mode == 'glm':

    print('Performing GLM..')
    glm.glm(args.dir, s=args.ws)
    print('')
    print('GLM comparisons carried out.')


else:
    error_msg()
















