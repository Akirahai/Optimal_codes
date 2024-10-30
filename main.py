import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data-dir', type=str, default='data/Knowledge_Base/', help='Directory for data dir')
    parser.add_argument('--q', type=int, default=2, help='Type of sequences') 
    parser.add_argument('--n', type=int, default=3, help='Length of sequences')
    parser.add_argument('--list-of-n',  type=int, nargs='+', default=[2,3,4,5,6], help='List of n')
    return parser.parse_args()
    
from libs import *



    
if __name__== "__main__":
    args = parse_args()
    q = args.q
    results = []
    # Make a tabulate for results comparison
    print("Starting calculations...")
    for n in tqdm(args.list_of_n, desc="Processing sequences", unit="sequence"):
        graph = QAryGraph(q, n)
        graph.build_graph()
        Optimal_codes, _ = graph.max_clique()
        
        if q == 2:
            vt_formula = VT_Formula(0, n)
            Assumed_codes = vt_formula.V_Ta_n()
        else:
            vt = VT(0, q, n, q*n)
            Assumed_codes = vt.find_solutions(0)
            
        results.append([q, n, Optimal_codes, Assumed_codes])
    
        # Print the results in tabular format
    headers = ["q", "n", "Optimal Codes", "Assumed Codes"]
    print(tabulate(results, headers=headers, tablefmt="pipe"))