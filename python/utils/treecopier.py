import hashlib
import multiprocessing
import os

import ROOT

from context_managers import open_root


ROOT.gROOT.SetBatch(True)

class TreeCopier(object):
    """
    Skim and copy the files of a sample according to a list of LFNs.

    Parameters
    ----------
    primary_dataset : str
                      The primary dataset name of the sample.
    LFN_path        : str
                      The path to the directory storing the LFN files.
    xrd_redirector  : str
                      The XRootD redirector prepended to the LFNs.
                      The default is the global redirector.
    cut             : str
                      The selection cut used to skim the sample.
    prep_out        : str
                      The local output directory storing the prepared samples.
    prefix          : str
                      The prefix for the sample subdirectory and merged file names.
    merge           : bool
                      Flag whether the files should be merged. The default is True.
    """
    def __init__(self, primary_dataset='', LFN_path='', xrd_redirector='cms-xrd-global.cern.ch', cut='', prep_out='', prefix='', merge=True):

        # Parse the file for the list of LFNs.
        LFN_file = os.path.join(LFN_path, primary_dataset + '.txt')
        with open(LFN_file, 'r') as infile:
            self.LFNs = [line.strip() for line in infile.readlines() if not line.isspace()]
        print 'Number of LFNs: {!s}\n'.format(len(self.LFNs))

        # Create the output directory for the copied files.
        self.copy_dir = os.path.join(prep_out, prefix + primary_dataset)
        try:
            os.makedirs(self.copy_dir)
            print 'Copy Destination: {}\n'.format(self.copy_dir)
        except OSError:
            if not os.path.isdir(self.copy_dir):
                raise

        # Set the remaining attributes.
        self.xrd_redirector = xrd_redirector
        self.cut = cut

        while self.LFNs:
            self._copy_files()

        if merge:
            self.merge_path = self.copy_dir + '.root'
            print 'Merged File: {}\n'.format(self.merge_path)
            self._merge_files()

    def _copy_files(self):
        """
        Copy the files in parallel, then remove all the
        successfully copied LFNs from the LFN list.
        """
        task_queue = multiprocessing.Queue()
        result_queue = multiprocessing.Queue()

        workers = [
            TreeCopyWorker(task_queue, result_queue, self.xrd_redirector, self.copy_dir, self.cut)
            for cpu in xrange(multiprocessing.cpu_count())
        ]

        for worker in workers:
            worker.start()

        for LFN in self.LFNs:
            task_queue.put(LFN)

        for worker in workers:
            task_queue.put(None)

        for worker in workers:
            worker.join()

        result_queue.put(None)
        for result in iter(result_queue.get, None):
            self.LFNs.remove(result)

    def _merge_files(self):
        """
        The first positional argument, isLocal, is set to false to prevent the
        files to merge from being copied to a temporary directory, as they are
        already available locally. PyROOT doesn't support keyword args.
        """
        merger = ROOT.TFileMerger(ROOT.kFALSE)
        merger.OutputFile(self.merge_path, 'RECREATE')

        for file_name in os.listdir(self.copy_dir):
            merger.AddFile(os.path.abspath(os.path.join(self.copy_dir, file_name)))

        merger.Merge()

class TreeCopyWorker(multiprocessing.Process):
    """
    A worker process which consumes LFNs from the task queue, returning
    them through the result queue if they were successfully copied.
    """
    def __init__(self, task_queue, result_queue, xrd_redirector, copy_dir, cut):
        super(TreeCopyWorker, self).__init__()
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.xrd_redirector = xrd_redirector
        self.copy_dir = copy_dir
        self.cut = cut

    def run(self):
        for LFN in iter(self.task_queue.get, None):
            input_url = 'root://{0}//{1}'.format(self.xrd_redirector, LFN)
            output_path = os.path.join(self.copy_dir, hashlib.sha256(LFN).hexdigest() + '.root')
            n_selected, n_total = self._copy_file(input_url, output_path)
            status = self._check_file(output_path)
            if status is not None:
                print '[PID {0!s}] Failure: {1}\n{2}\n'.format(os.getpid(), status, LFN)
            else:
                print '[PID {0!s}] Success: Selected {1!s} out of {2!s} entries.\n{3}\n'.format(os.getpid(), n_selected, n_total, LFN)
                self.result_queue.put(LFN)

    def _copy_file(self, input_url, output_path):
        """
        Copy the file at the remote url, saving a skimmed copy locally.
        Objects other than the TTree 'tree' are copied directly.
        """
        with open_root(input_url, 'r') as infile, open_root(output_path, 'w'):
            for key in infile.GetListOfKeys():
                if key.GetName() == 'tree':
                    continue
                obj = key.ReadObj()
                obj.Write()
            intree = infile.Get('tree')
            n_total = intree.GetEntriesFast()
            outtree = intree.CopyTree(self.cut)
            n_selected = outtree.GetEntriesFast()
            outtree.Write()
        return n_selected, n_total

    def _check_file(self, output_path):
        """
        Check that the file is not corrupted.
        """
        try:
            with open_root(output_path, 'r') as outfile:
                if not outfile.GetNkeys():
                    return 'Missing TKeys'
                elif outfile.TestBit(ROOT.TFile.kRecovered):
                    return 'kRecovered'
                elif outfile.IsZombie():
                    return 'kZombie'
                else:
                    return None
        except IOError:
            return 'Unable to open for error checking.'

