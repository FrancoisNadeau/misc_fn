#!/usr/bin/python3

import json
import os
import regex
import bibtexparser
import bs4
import re

from argparse import ArgumentParser
from bs4 import BeautifulSoup
from os import PathLike
from pathlib import Path, PosixPath
from typing import Union
from unidecode import unidecode

from get_desc import get_desc

def stripstr(txtlist: list) -> str:
    """
    Returns whitespace-stripped list of strings found in ``txtlist``.
    
    Only the leading or trailing whitespace characters are removed.
    """

    return list(map(str.strip, txtlist))


def rep2dict(pattern: Union[str, re.Pattern],
             txt: str
             ) -> dict:
    
    """
    Returns a dict mapped as {strings found: text inbetween}.
    
    Args:
        pattern: str or re.Pattern
            Pattern to look for in text.
            Must be compliant with the ``re`` python3 module.
        txt: str
            String representation of the text to be divided into sections.
    
    Returns: dict
    """

    keys = [itm.replace(':', '').strip() for itm in
            re.findall(pattern, txt)]
    sections = tuple(re.split(pattern, txt))[1:]
    keys, sections = stripstr(keys), stripstr(sections)
    return dict(zip(keys, sections))


def parse_func(report: dict,
               func_pat: Union[str, re.Pattern],
               subtitle_pat: Union[str, re.Pattern],
               fsub_pat: Union[str, re.Pattern],
               short: bool = True
               ) -> dict:
    """
    Returns FMRIPrep report "Functional" section as dict.
    """

    # Functional section
    test = report['Functional']
    # Functional report headers
    fkeys = re.findall(func_pat, test)
    # Functional report subtitles and text within
    fvalues = rep2dict(func_pat, test)
    ckey, skey = 'Confounds collected', 'Summary'
    [fvalues.update({s[0]: rep2dict(fsub_pat, s[1])})
                   for s in tuple(fvalues.items())]
    trans = {'(': '', ')': '', ' ': '_'}
    [fvalues[key].update({**{[itm[0].translate(str.maketrans(trans))+'_s'
                              if itm[1].endswith('s') else
                              itm[0].translate(str.maketrans(trans))][0]:
                             re.sub(':', '',
                                    itm[1]).split('\n'*3)[0].rstrip('s')
                             for itm in
                             tuple(rep2dict(subtitle_pat,
                                     fvalues[key][skey]).items())},
                          **{ckey: fvalues[key][ckey].split('\n'*4)[0][:-1]}})
     for key in tuple(fvalues.keys())]
    [fvalues[key].pop(skey) for key in tuple(fvalues.keys())]
    if short is True:
        [fvalues[key].pop(ckey) for key in tuple(fvalues.keys())]
    return fvalues


def parse_summary(report: dict,
                  subtitle_pat: Union[str, re.Pattern]
                  ) -> dict:
    """
    Returns FMRIPrep report "Summary" section as dict.
    """
    
    test = report['Summary'].replace('Task: ', 'Task-')
    fkey = 'Functional series'
    test = rep2dict(subtitle_pat, test)
    test.update({'Tasks': test[fkey].split('\n\n')[1].split('\n')})
    test.update({fkey: test[fkey].split('\n\n')[0]})
    test['Subject ID'] = 'sub-'+test['Subject ID']
    sos = test['Standard output spaces']
    test['Standard output spaces'] = list(map(str.strip, sos.split(',')))
    nsos = test['Non-standard output spaces']
    test['Non-standard output spaces'] = list(map(str.strip, nsos.split(',')))
    return test


def parse_boilerplate(report: dict,
                      bp_pat: Union[str, re.Pattern]
                      ) -> dict:
    """
    Returns FMRIPrep report "boilerplate" section attributes as dict.
    """

    test = report['boilerplate'].split('\n\n', 1)[1]
    test = rep2dict(bp_pat, test)
    [report.update({itm[0]: re.sub('\n{2,}#*', '',
                                   itm[1].replace(': ', '').strip())})
     for itm in tuple(test.items())]
    report.pop('boilerplate'), report.pop('References')
    return report


def parse_bib(report: dict) -> dict:
    """
    Returns a dict containing entries from the "Bibliography" section.
    """

    from bibtexparser import loads as bibloads

    return dict(enumerate(bibloads(report['Bibliography']).entries))


def make_ses_ids(report: dict) -> list:
    """
    Returns a list of BIDS compliant identifiers for each fMRI session.
    """

    base = tuple(report['Functional'].keys())
    sessions = ['_'.join((report['Summary']['Subject ID'],
                          'ses-'+b.split('session ')[1].split(',')[0],
                          'task-'+b.split('task ')[1][:-1]))
                for b in base]
    return sessions


def parse_report(src: Union[str, PathLike, PosixPath],
                 json_file: Union[str, PathLike, PosixPath] = None,
                 # title_pat: Union[str, re.Pattern] = title_pat,
                 # func_pat: Union[str, re.Pattern] = func_pat,
                 # subtitle_pat: Union[str, re.Pattern] = subtitle_pat,
                 features: str = 'lxml',
                 ensure_ascii: bool = True,
                 **kwargs
                 ) -> dict:
    """
    Returns a dict representation of an FMRIPrep report.
    
    Args:
        src: str, PathLike or PosixPath
            Path to the html report file.

        json_file: str, PathLike or PosixPath (Default = None)
            Path to a .json file containing user-defined regex
            pattern alternatives to parse the html report.

        features: str (Default = 'lxml')
            String representing the desired html parser engine.
            See ``help(bs4.BeautifulSoup.__init__)`` for details.

        ensure_ascii: bool (Default = True)
            Boolean value indicating if non-ASCII characters
            should be converted accordingly.
            Flawless conversion is done ``unidecode.unidecode``,
            which supports multiple languages with various alphabet types.
    
    Returns: dict
        Dict representation of an FMRIPrep report.
        The outermost keys are
        ['Summary', 'Anatomical', 'Functional', 'About', 'problems',
         'figures', 'Bibliography', 'Errors', 'Anatomical_preprocess',
         'Functional_preprocess', 'Copyright'].
    """

    from bs4 import BeautifulSoup as bsf
    from unidecode import unidecode
    from os.path import isfile

# Compile regex patterns
    progdir = os.path.dirname(__file__)
    if json_file is None:
        json_file = os.path.join(progdir, 'ReportParserPatterns.json')
    with open(json_file, 'r') as jfile:
        patterns = json.load(jfile)
        jfile.close()

    err_pat, getfig_pat, func_pat, fsub_pat, \
    title_pat, subtitle_pat, bp_pat = \
    tuple(map(re.compile, patterns.values()))

    if os.path.isfile(src):
        htmltext = Path(src).read_text().strip()
    htmltext = htmltext[re.search(title_pat, htmltext).span()[0]:]
    keys = [re.search('\"(\w*?)\"', t).group(1)
            for t in re.findall(title_pat, htmltext)]
    txt_base = enumerate(list(re.split(title_pat, htmltext))[1:])
    txt_values = [bsf(t[1], features).text.strip(keys[t[0]]).strip()
                  for t in txt_base]
    figures = [[fig.strip() for fig in getfig_pat.findall(t)]
               for t in txt_values]
    error_messages = [[err.strip() for err in err_pat.findall(t)]
                      for t in txt_values]
    txt_values = [err_pat.sub('', t) for t in txt_values]
    txt_values = [getfig_pat.sub('', err_pat.sub('', t))
                  for t in txt_values]
    text_ = dict(zip(keys, txt_values))
    errs_ = dict(zip(keys, error_messages))
    figs_ = dict(zip(keys, figures))
    parsed_ = text_
    parsed_['Problems'] = errs_
    parsed_['Figures'] = figs_
    # Polishing report's primary sections
    parsed_['Summary'] = parse_summary(parsed_, subtitle_pat)
    parsed_['Anatomical'] = rep2dict('.*\: ', parsed_['Anatomical'])
    n_discarded = parsed_['Anatomical']['Discarded images'].split('\n')[0]
    parsed_['Anatomical']['Discarded images'] = n_discarded
    parsed_['Functional'] = parse_func(parsed_, func_pat,
                                       subtitle_pat, fsub_pat)
    sessions = make_ses_ids(parsed_)
    parsed_['Functional'] = {sessions[val[0]]: val[1] for val in
                             enumerate(tuple(parsed_['Functional'].values()))}
    parsed_ = parse_boilerplate(parsed_, bp_pat)
    parsed_['About'] = rep2dict(subtitle_pat, parsed_['About'])
    parsed_['Bibliography'] = parse_bib(parsed_)
    parsed_['Errors'] = parsed_['errors'].splitlines()[1:]
    parsed_.pop('errors')
    # Setting keys without whitespaces
    adp = parsed_['Anatomical data preprocessing']
    parsed_['Anatomical_preprocess'] = adp.strip()
    parsed_.pop('Anatomical data preprocessing')
    fdp = parsed_['Functional data preprocessing']
    parsed_['Functional_preprocess'] = fdp.replace('\n\n\n###', '').strip()
    parsed_.pop('Functional data preprocessing')
    parsed_['Copyright'] = parsed_['Copyright Waiver']
    parsed_.pop('Copyright Waiver')
    method_keys = ['Anatomical_preprocess', 'Functional_preprocess',
                   'Errors', 'Copyright']
    # Extracting metrics units
    ovs = parsed_['Anatomical']['Output voxel size']
    unit = re.search(r'\d{1}(.*?)\s{1}', ovs).group(1)
    ovs = re.sub(unit, '', str(tuple(ovs.split(f'{unit} x '))))
    parsed_['Anatomical'].update({f'Output voxel size {unit}': ovs})
    parsed_['Anatomical'].pop('Output voxel size')
    od = str(tuple(parsed_['Anatomical']['Output dimensions'].split('x')))
    parsed_['Anatomical'].update({'Output dimensions': od})
    parsed_['Methods'] = {key: parsed_[key] for key in method_keys}
    [parsed_.pop(key) for key in method_keys]
    [parsed_.update({key: {itm[0].replace(' ', '_'): itm[1]
                           for itm in tuple(parsed_[key].items())}})
     for key in tuple(parsed_.keys())[:-2]]

    return parsed_

def main():
    description, help_msgs = get_desc(parse_report)
    parser_ = ArgumentParser(prog=parse_report,
                             desc=description.splitlines()[0],
                             usage=description)
    parser_.add_argument('src', nargs=1, help=help_msgs[0])
    parser_.add_argument(*('-j', '--json-file'), dest='json_file',
                        nargs='*', help=help_msgs[1])
    parser_.add_argument(*('-f', '--features'), dest='features', default='lxml',
                        nargs='*', help=help_msgs[2])
    parser_.add_argument('--ensure-ascii', dest='ensure_ascii', type=bool,
                        action='store_true', help=help_msgs[-1])
    args = parser_.parse_arguments()
    parse_report(args.src[0], args.json_file,
                 args.features, args.ensure_ascii)

if __name__ == '__main__':
    main()
