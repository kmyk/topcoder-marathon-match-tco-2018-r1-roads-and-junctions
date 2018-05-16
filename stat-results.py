# Python Version: 3.x
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import collections
import json
from tabulate import tabulate  # https://pypi.org/project/tabulate/

def load_list_of_json_file(path):
    df = []
    with open(path) as fh:
        decoder = json.JSONDecoder(object_pairs_hook=collections.OrderedDict)
        for i, line in enumerate(fh):
            df += [ decoder.decode(line) ]
    df = pd.DataFrame(df, columns=df[0].keys())
    df.set_index('seed', inplace=True)
    return df

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('what', choices='table summary compare pairplot distplot')
    parser.add_argument('file', nargs='?', default='/dev/stdin')
    parser.add_argument('--seed', type=int)
    parser.add_argument('--compare')
    parser.add_argument('--save')
    args = parser.parse_args()

    df = load_list_of_json_file(args.file)

    if args.what == 'table':
        headers = [ df.index.name ] + list(df.columns)
        for key in list(headers):
            if key.endswith('samples'):
                df = df.drop(key, axis=1)
                headers.remove(key)
        df = df.sort_index()
        s = tabulate(df, headers=headers, showindex='always', tablefmt='orgtbl')
        lines = s.splitlines()
        lines[1] = lines[1].replace('+', '|')
        s = '\n'.join(lines)
        print(s)

    elif args.what == 'summary':
        print('average of average reference delta =', df['average_reference_delta'].mean())
        print('median of average reference delta =', df['average_reference_delta'].median())
        print('minimum of average reference delta =', df['average_reference_delta'].min())

    elif args.what == 'compare':
        if args.compare is None:
            parser.error('the following arguments are required: --compare')
        df1 = df
        df2 = load_list_of_json_file(args.compare)
        df1 = df1.rename(columns={ 'NJ': 'NJ1', 'average_reference_delta': 'ave_delta_1' })
        df2 = df2.rename(columns={ 'NJ': 'NJ2', 'average_reference_delta': 'ave_delta_2' })
        df = df1.join(df2[ [ 'NJ2', 'ave_delta_2' ] ])
        df = df.assign(ave_delta_diff=lambda row: row.ave_delta_1 - row.ave_delta_2)
        df = df.sort_values(by='ave_delta_diff')
        headers = [ 'seed', 'S', 'NC', 'junction_cost', 'failure_probability', 'reference_score', 'NJ1', 'ave_delta_1', 'NJ2', 'ave_delta_2', 'ave_delta_diff' ]
        for key in list(df.columns):
            if key not in headers:
                df = df.drop(key, axis=1)
        s = tabulate(df, headers=headers, showindex='always', tablefmt='orgtbl')
        lines = s.splitlines()
        lines[1] = lines[1].replace('+', '|')
        s = '\n'.join(lines)
        print(s)

    elif args.what.endswith('plot'):
        if args.what == 'pairplot':
            sns.pairplot(df, vars='S NC junction_cost failure_probability reference_score NJ average_reference_delta'.split())
        elif args.what == 'distplot':
            if args.seed is not None:
                try:
                    sns.distplot(df.ix[args.seed]['delta_samples'])
                except np.linalg.linalg.LinAlgError:
                    sns.distplot(df.ix[args.seed]['delta_samples'], kde=False)
            else:
                sns.distplot(df['average_reference_delta'], norm_hist=False, vertical=True, kde=False)
        else:
            assert False
        if args.save:
            plt.savefig(args.save)
        else:
            plt.show()

    else:
        assert False


if __name__ == '__main__':
    main()
