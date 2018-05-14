# Python Version: 3.x
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import collections
import json
from tabulate import tabulate  # https://pypi.org/project/tabulate/

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('what', choices='table summary pairplot distplot')
    parser.add_argument('file', nargs='?', default='/dev/stdin')
    parser.add_argument('--seed', type=int)
    parser.add_argument('--save')
    args = parser.parse_args()

    # load data from list-of-json format
    df = []
    with open(args.file) as fh:
        decoder = json.JSONDecoder(object_pairs_hook=collections.OrderedDict)
        for i, line in enumerate(fh):
            df += [ decoder.decode(line) ]
    df = pd.DataFrame(df, columns=df[0].keys())
    df.set_index('seed', inplace=True)

    # plot
    if args.what == 'table':
        headers = list(df.columns)
        for key in list(headers):
            if key.endswith('samples'):
                df = df.drop(key, axis=1)
                headers.remove(key)
        s = tabulate(df, headers=headers, showindex='always', tablefmt='orgtbl')
        lines = s.splitlines()
        lines[1] = lines[1].replace('+', '|')
        s = '\n'.join(lines)
        print(s)
    elif args.what == 'summary':
        print('average of average reference delta =', df['average_reference_delta'].mean())
        print('median of average reference delta =', df['average_reference_delta'].median())
        print('minimum of average reference delta =', df['average_reference_delta'].min())
    elif args.what.endswith('plot'):
        if args.what == 'pairplot':
            sns.pairplot(df, vars='S NC junction_cost failure_probability reference_score NJ average_reference_delta'.split())
        elif args.what == 'distplot':
            if args.seed is None:
                parser.error('the following arguments are required: --seed')
            try:
                sns.distplot(df.ix[args.seed]['delta_samples'])
            except np.linalg.linalg.LinAlgError:
                sns.distplot(df.ix[args.seed]['delta_samples'], kde=False)
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
