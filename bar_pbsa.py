import yaml
import numpy as np
import pandas as pd
import argparse
import logging
import parmed.amber
import parmed.tools
from pathlib import Path
import pytraj as pt
import shutil


def init_logger(log_path='mbar_pbsa.log'):

    logging.basicConfig(level=logging.INFO,
            format='%(asctime)s : %(levelname)s : %(message)s',
            filename=log_path)
    logging.getLogger().addHandler(logging.StreamHandler())


def parse_args():

    '''
    cli to 1) strip trajectories
           2) prepare BAR/PBSA input parameters
           3) run BAR/PBSA in parallel
           4) post-process BAR data
    '''

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    parser_a = subparsers.add_parser('strip')
    parser_a.add_argument('strip_input', type=str, help='path to input strip yaml file')

    parser_b = subparsers.add_parser('prep')
    parser_b.add_argument('prep_input', type=str, help='path in input prep yaml file')

    return parser.parse_args()


class strip_traj:

    def __init__(self, yaml_path):

        self.yaml_path = yaml_path
        self.parse_input()


    def parse_input(self):

        '''parse input yaml'''

        with open(self.yaml_path) as f:
            inp = yaml.safe_load(f)
            logging.info(inp)

        self.dest_path = inp['dest_path']
        self.complex_paths = inp['complex_paths']
        self.ligand_paths = inp['ligand_paths']
        self.ligand_mask = inp['ligand_mask']
        self.ion_decharge = inp['ion_decharge']
        self.last_half_frames = inp['last_half_frames']


    def get_decharged_ion(self, run_path):

        '''find index of decharged ion to keep'''

        assert self.ion_decharge

        # load last parm at 1.000
        last_parm = list(Path(run_path, '1.000').glob('*.parm7'))[0]
        parm_in = parmed.amber.AmberParm(str(last_parm))

        # find Na+ or Cl- with 0 charge
        ions = parmed.tools.printDetails(parm_in, ':Cl-,:Na+')
        ions.execute()
        ions_str = str(ions)

        for line in ions_str.split('\n'):
            if len(line.split()) > 9 and \
               line.split()[9] not in ['-1.0000', '1.0000'] and \
               'ATOM' not in line:

                ion_line = line
                break

        # get index
        atom_idx = ion_line.split()[0]

        return atom_idx


    def make_folds(self):

        '''make folders to save stripped trajectories'''

        self.dest_dir = Path(self.dest_path)

        # make for stripped complex and ligands
        self.com_folds = [f't{i}' for i in range(1, len(self.complex_paths) + 1)]
        self.lig_folds = [f't{i}' for i in range(1, len(self.ligand_paths) + 1)]

        self.lig_lamdas = sorted([x.name for x in Path(self.ligand_paths[0]).iterdir()])
        self.com_lamdas = sorted([x.name for x in Path(self.complex_paths[0]).iterdir()])

        for com_fold in self.com_folds:
            for com_lamda in self.com_lamdas:
                Path(self.dest_dir, 'complex', com_fold, com_lamda).mkdir(parents=True, exist_ok=True)
                logging.info(f"mkdir {Path(self.dest_dir, 'complex', com_fold, com_lamda)}")

        for lig_fold in self.lig_folds:
            for lig_lamda in self.lig_lamdas:
                Path(self.dest_dir, 'ligands', lig_fold, lig_lamda).mkdir(parents=True, exist_ok=True)
                logging.info(f"mkdir {Path(self.dest_dir, 'ligands', lig_fold, lig_lamda)}")

        # for concat traj
        for com_lamda in self.com_lamdas:
            Path(self.dest_dir, 'concat_complex', com_lamda).mkdir(parents=True, exist_ok=True)
            logging.info(f"mkdir {Path(self.dest_dir, 'concat_complex', com_lamda)}")

        for lig_lamda in self.lig_lamdas:
            Path(self.dest_dir, 'concat_ligands', lig_lamda).mkdir(parents=True, exist_ok=True)
            logging.info(f"mkdir {Path(self.dest_dir, 'concat_ligands', com_lamda)}")


    def strip_single_traj(self, ligcom, run_path, lamda):

        '''strip waters and ions, use amber mask to select which atoms to keep'''

        assert ligcom in ['ligands', 'complex']

        parm = str(list(Path(run_path, lamda).glob('*.parm7'))[0])
        trajin = str(list(Path(run_path, lamda).glob('*.nc'))[0])

        traj_name = trajin.split('/')[-1]
        parm_name = parm.split('/')[-1]

        traj = pt.load(trajin, parm)
        n_frames = traj.n_frames

        # strip 
        if ligcom == 'ligands':
            # ligand only
            mask = f'{self.ligand_mask}'
            paths = self.ligand_paths

        elif ligcom == 'complex':
            # ligand and protein
            # residues from number of ca
            n_residues = len(traj.top.select('@CA'))
            mask = f':1-{n_residues}|{self.ligand_mask}'
            paths = self.complex_paths
            
        # get ion index
        if self.ion_decharge:
            atom_idx = self.get_decharged_ion(run_path)
            mask = mask + f'|@{atom_idx}'

        # select region to keep
        traj = traj[mask]
        logging.info(f'{run_path}, lambda: {lamda}, mask: {mask}')

        # last half
        if self.last_half_frames:
            traj = traj[n_frames//2:]

        # save nc and restart
        path_idx = paths.index(run_path)

        if ligcom == 'ligands':
            fold = self.lig_folds[path_idx]
            out_traj = str(Path(self.dest_dir, 'ligands', fold, lamda, traj_name))
            out_parm = str(Path(self.dest_dir, 'ligands', fold, lamda, parm_name))

        elif ligcom == 'complex':
            fold = self.com_folds[path_idx]
            out_traj = str(Path(self.dest_dir, 'complex', fold, lamda, traj_name))
            out_parm = str(Path(self.dest_dir, 'complex', fold, lamda, parm_name))

        pt.write_traj(out_traj, traj)
        pt.write_parm(out_parm, traj.top)


    def align(self, ligcom, lamda):

        '''concat all replicate trajectories, rmsd fit to first frame and autoimage'''

        # glob all trajs
        if ligcom == 'ligands':
            paths = sorted(list(Path(self.dest_dir, 'ligands').glob(f'*/{lamda}/*.nc')))
            parm = sorted(list(Path(self.dest_dir, 'ligands').rglob('*.parm7')))[0]
        else:
            paths = sorted(list(Path(self.dest_dir, 'complex').glob(f'*/{lamda}/*.nc')))
            parm = sorted(list(Path(self.dest_dir, 'complex').rglob('*.parm7')))[0]

        paths = list(map(str, paths))

        # read first parm
        parm = str(parm)

        # iterload all
        traj = pt.iterload(paths, parm)

        # rmsd align to first
        if ligcom == 'ligands':
            mask = self.ligand_mask
        elif ligcom == 'complex':
            n_residues = len(traj.top.select('@CA'))
            mask = f':1-{n_residues}@CA'

        rmsd = pt.rmsd(traj, mask=mask, ref=0)

        # autoimage
        traj = traj.autoimage()

        # save concat traj, restart, parm
        # bug? only works with ncrst restarts, not rst7
        traj_save = str(Path(self.dest_dir, f'concat_{ligcom}', lamda, 'ti001.nc'))
        restart_save = Path(self.dest_dir, f'concat_{ligcom}', lamda, 'ti.ncrst')
        parm_save = str(Path(self.dest_dir, f'concat_{ligcom}', lamda, 'ti.parm'))

        pt.write_traj(traj_save, traj)
        pt.write_traj(str(restart_save), traj, frame_indices=[0])
        pt.save(parm_save, traj.top)

        # fix restart name, remove suffix
        Path(f'{restart_save}.1').rename(restart_save)

        logging.info(f'align {self.dest_path} {ligcom} {lamda}')


    def strip_all_traj(self):

        for lig_path in self.ligand_paths:
            for lamda in self.lig_lamdas:
                self.strip_single_traj('ligands', lig_path, lamda)

        for com_path in self.complex_paths:
            for lamda in self.com_lamdas:
                self.strip_single_traj('complex', com_path, lamda)


    def align_all_traj(self):

        for lig_lambda in self.lig_lamdas:
            self.align('ligands', lig_lambda)

        for com_lambda in self.com_lamdas:
            self.align('complex', com_lambda)


    def run_all(self):

        self.make_folds()
        self.strip_all_traj()
        self.align_all_traj()


class prep_bar_pbsa:

    def __init__(self, yaml_path):

        self.yaml_path = yaml_path
        self.parse_input()


    def parse_input(self):

        '''parse input yaml'''

        with open(self.yaml_path) as f:
            inp = yaml.safe_load(f)
            logging.info(inp)

        self.dest_path = inp['dest_path']
        self.ligand_res = inp['ligand_res']
        self.istrng = inp['istrng']
        self.epsin = inp['epsin']
        self.radiscale = inp['radiscale']
        self.protscale = inp['protscale']


    def get_n_frames(self, ligcom):

        '''frame counts from first traj'''

        assert ligcom in ['ligands', 'complex']

        trajin = sorted(list(Path(self.dest_path, f'concat_{ligcom}').rglob('*.nc')))[0]
        parm = sorted(list(Path(self.dest_path, f'concat_{ligcom}').rglob('*.parm')))[0]

        traj = pt.iterload(str(trajin), str(parm))

        return traj.n_frames


    def write_sander(self, ligcom):

        # set fillratio based on ligcom
        if ligcom == 'ligands':
            fillratio = 4
            dest = Path(self.lig_dest, 'pb_input.txt')
        elif ligcom == 'complex':
            fillratio = 2
            dest = Path(self.com_dest, 'pb_input.txt')

        n_frames = self.get_n_frames(ligcom)

        template = f'''PB energy calculation
 &cntrl
   imin = 5, nstlim = {n_frames}, dt = 0.002,
   ntx = 2, irest = 0,
   ntb = 0,
   ntt = 1,
   temp0 = 298.0, ig = -1,
   ntp = 0,
   ntc = 2, ntf = 2,
   ioutfm = 1, iwrap = 0,
   ntwe = 1000, ntwx = 10000, ntpr = 1, ntwr = 10000,
   cut = 9999.0,
   maxcyc = 1,
   ipb=2, inp=0

/
&pb
   npbverb=0, istrng={self.istrng}, epsout=80.0, epsin={self.epsin}, space=.5,
   accept=0.001, radiopt=0, fillratio={fillratio},
   npbopt=0, bcopt=1, solvopt=1, maxitn=10000,
   frcopt=0, nfocus=2, radiscale={self.radiscale}, radires="{self.ligand_res}", protscale={self.protscale},
   dprob=1.4
/
'''

        Path(dest).write_text(template)
        logging.info(f'{self.dest_path} {ligcom}\n{template}')


    def make_folds(self):

        self.lig_path = Path(self.dest_path, 'param_sweep_ligands')
        self.com_path = Path(self.dest_path, 'param_sweep_complex')

        folder = f'e_{self.epsin}_r_{self.radiscale}_p_{self.protscale}'

        self.lig_dest = Path(self.lig_path, folder)
        self.com_dest = Path(self.com_path, folder)

        # copy clean dir over
        com_source = Path(self.dest_path, 'concat_complex')
        lig_source = Path(self.dest_path, 'concat_ligands')

        shutil.copytree(com_source, str(self.com_dest))
        shutil.copytree(lig_source, str(self.lig_dest))

        logging.info(f'setup {self.lig_dest}')
        logging.info(f'setup {self.com_dest}')


    def run_all(self):

        self.make_folds()
        self.write_sander('ligands')
        self.write_sander('complex')


# write pbsa input files, set epsin, radiscale, protscale

# multiprocess sander runs

# bar calculations


if __name__ == '__main__':

    init_logger()
    args = parse_args()

    #if args.input:
    #    test = strip_traj(args.input)
    #    test.run_all()

    if args.prep_input:
        test = prep_bar_pbsa(args.prep_input)
        test.run_all()
