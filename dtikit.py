import os
import logging
import shelve
import re
import gzip
from io import StringIO
from collections import namedtuple
import pandas as pd
import requests
from bioservices.uniprot import UniProt
from Bio import SwissProt

class Chemical():
    __slots__ = ('cid',
    'stereo',
    'name',
    'interactions',
    'inhibits',
    'activates',
    'molecular_weight',
    'SMILES_string')

    def __init__(self,cid,**kwargs):
        self.cid = cid
        for k,v in kwargs.iteritems():
            setattr(self, k,v)

    def __repr__(self):
        return '<Chemical: {0}>'.format((self.name if hasattr(self,'name') else self.cid))


class Protein():
    __slots__ =[ 'features',
    'sequence',
    'cross_references',
    'host_organism',
    'references',
    'molecule_type',
    'keywords',
    'host_taxonomy_id',
    'organism_classification',
    'protein_existence',
    'comments',
    'seqinfo',
    'description',
    'taxonomy_id',
    'sequence_length',
    'created',
    'annotation_update',
    'organelle',
    'accessions',
    'entry_name',
    'sequence_update',
    'organism',
    'gene_name',
    'data_class',
    'interactions',
    'inhibitors',
    'activators',
    'names']

    def __init__(self,entry_name,**kwargs):
        self.entry_name = entry_name
        for k,v in kwargs.iteritems():
            setattr(self,k,v)

    def __repr__(self):
        return '<Protein: {0}>'.format(self.entry_name)


class DTIkit():
    "Drug-target interaction data toolkit"

    def __init__(self, folder=None, species='9606', ignore_stereo=False, stitch_version='5.0', db_version='5.0.7',\
            stitch_channel='combined_score', rebuild=False):
        """Create a DTIkit object that enables easy access to drug, target, and interaction data.

        Parameters
        ----------
        *folder* : str, default None
            Folder with the DrugBank and STITCH files. If not specified, and if this package was used before, 
            remembers and uses the last folder input (i.e. the last used location).
        *species* : str, default "9606"
            the NCBI taxonomy ID of the species (ex: 9606 for human)
        *ignore_stereo* : boolean, default False
            tells the parser to treat stereoisomers of a chemical as equivalent
        *stitch_version* : string, default "5.0"
            version of the STITCH database to use
        *db_version* : string, default "5.0.7"
            version of the DrugBank database to use
        *stitch_channel* : string, default "combined_score"
            the stitch information channel (e.g. 'experimental' or 'database' etc.)
        *rebuild* : boolean, default False
            whether to rebuild the internal data structures even if they already exist
        """
        if folder is None:
            if _config['last_folder'] is None:
                raise ValueError('Package requires providing the folder with DrugBank and STITCH raw files at least once')
            else:
                folder = _conf['last_folder']
        else:
            _config['last_folder'] = folder
        self.folder = folder
        if species == 'all':
            species=''
        else:
            species=species+'.'
        suffix = '_'.join([species, 'stereo'+str(ignore_stereo), 
            'STv'+stitch_version+'('+stitch_channel+')', 'DBv'+db_version])+'.shl'
        self._files = dict([(element,folder+element+suffix) for element in ['proteins','chemicals','names']])
        for f in self._files.values():
            if not os.path.isfile(f):
                rebuild=True
        for k,v in self._files.iteritems():
            setattr(self,k,shelve.open(filename=v, flag=('w' if rebuild else 'r'), protocol=2, writeback=True)
        if rebuild:
            self._build(species=species, ignore_stereo=ignore_stereo, stitch_version=stitch_version,
                db_version=db_version, stitch_channel=stitch_channel)


    def _build(self,species,ignore_stereo,stitch_version,db_version,stitch_channel):
        "Builds the internal data structures"
        logging.warning('Building data structures. This is done only once but takes quite a while. Please be patient.')
        links_file = self.folder+os.sep+species+"protein_chemical.links.detailed.v"+stitch_version+".tsv.gz"
        if not os.path.isfile(links_file):
            links_file = pd.folder+os.sep+species+"protein_chemical.links.v"+stitch_version+".tsv.gz"
            if stitch_channel is not 'combined_score':
                raise IOError('stitch_channel is {0} but detailed links file not found:\n{1}'.format(stitch_channel,links_file))
        logging.info('Reading STITCH files')
        links = pd.read_csv(links_file,'\t')
        u = UniProt(verbose=False)
        #chems = pd.read_csv(self.folder+os.sep+'chemicals.v'+stitch_version+'.tsv.gz','\t').set_index('chemical')
        logging.info('Processing each protein-chemical link')
        for row in links.itertuples():
            if row.protein not in self.proteins:
                self._proc_protein(row.protein)
            if row.chemical not in self.chemicals:
                self.chemicals[row.chemical] = Chemical(row.chemical)
            self.proteins[row.protein].interactions.append(row.chemical)
            self.chemicals[row.chemical].interactions.append(row.protein)
#        self._proc_all_chemicals()


#    def _proc_all_chemicals(self):
#        logging.info('Processing all chemicals for structure')
#        files = [self.folder+os.sep+'chemicals.v'+stitch_version+'.tsv.gz',
#
#
#        with gzip.open(self.folder+os.sep+'chemicals.v'+stitch_version+'.tsv.gz') as f:
#            headers = f.readline().strip().split('\t')
#            for line in f:
#                row = dict(zip(headers,line.strip().split('\t')))
#                if row['chemical'] in self.chemicals:
#                    for k,v in row.iteritems():
#                        if k is not 'chemical':
#                            setattr(self.chemicals[row['chemical']],k,v)
#                    self.names[row['name']]=row['chemical']
#        with gzip.open(self.folder+os.sep+'chemical.aliases.'+stitch_version+'.tsv.gz') as f:


    def _proc_protein(self,protid):
        logging.info('Processing '+protid)
        org,prot = protid.split('.')
        s = u.search('taxonomy:'+org+'+and+'+prot,frmt='txt')
        if len(s) > 0:
            prots = SwissProt.parse(StringIO(s))
            prot = prots[0]
            for p in prots: # If other matches, pick a reviewed one if possible
                if p.data_class=='Reviewed':
                    prot = p
                    break
            p = Protein(names=[],interactions=[],inhibitors=[],activators=[],**prot.__dict__)
            exp = re.compile('( )?(RecName: |AltName: )?(Full=|Short=){1}(?P<name>[^\{\}]+)')
            for s in p.description[:-1].split(';'):
                m = exp.match(s)
                if m is not None:
                    name = m.group('name')
                    p.names.append(name)
                    self.names[name] = protid
            self._proteins[protid] = p # Record the protein
            for acc in p.accessions:   # Now record the names of the protein
                self.names[acc] = protid 
            self.names[p.entry_name] = protid
            
