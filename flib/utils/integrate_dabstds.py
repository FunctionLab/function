from optparse import OptionParser
from collections import defaultdict
import random
import os
import sys

# CONSTANTS
binaries = {'Counter': 'Counter',
            'NetworkCombiner': 'NetworkCombiner',
            'DChecker': 'DChecker',
            'Dat2Dab': 'Dat2Dab'}

# FUNCTIONS
# functions that call counter (3 stages)
# Counter Learn


def learn(job=None, job_name=None, global_answers=None, holdout=None,
          counter=None, data=None, working=None, zeros=None, standards=None,
          stddir=None, wtsdir=None, weights=None, flipneg=None,
          multiplier=None, noweightneg=None, local_cmds=None):

    # Make global network.
    learn_jobs = {}
    try:
        os.mkdir(working + '/counts')
    except OSError:
        pass

    cmdline = counter + ' -w ' + global_answers + \
        ' -d ' + data + \
        ' -o ' + working + '/' + \
        ' -Z ' + zeros

    if job:
        job.set_name_command(job_name + '-global-learn', cmdline)
        learn_jobs.append(job.submit(os.path.join(working, 'global-learn.sh')))
    elif local_cmds is not None:
        local_cmds['global'] += cmdline + ';'

    for standard in standards:
        if os.path.exists(stddir + '/' + standard + '.dab'):
            cmdline = counter + ' -w ' + stddir + '/' + standard + '.dab' + \
                ' -d ' + data + \
                ' -o ' + working + '/counts/' + \
                ' -O ' + standard + \
                ' -Z ' + zeros
        else:
            cmdline = counter + ' -w ' + stddir + '/' + standard + '.dat' + \
                ' -d ' + data + \
                ' -o ' + working + '/counts/' + \
                ' -O ' + standard + \
                ' -Z ' + zeros
        if wtsdir is not None:
            cmdline += ' -W ' + weights[standard]
            if flipneg is not True:
                cmdline += ' -F'
            if multiplier is not None:
                cmdline += ' -f ' + str(multiplier)
            if noweightneg is True:
                cmdline = cmdline + ' -N'

        if holdout is not None:
            cmdline += ' -G ' + holdout

        if job:
            job.set_name_command(job_name + '-' + standard + '-learn', cmdline)
            learn_jobs.append(job.submit(os.path.join(working, standard
                                                      + '-learn.sh')))
        elif local_cmds is not None:
            local_cmds[standard] += cmdline + ';'

    if job:
        job.set_depends(None)
    return learn_jobs

# Counter Networks


def networks(job=None, job_name=None, counter=None, data=None, working=None,
             alphas=None, pseudo=None, stddir=None, depends=None,
             weights_suffix=None):
    if depends is not None:
        job.set_depends(depends[:])
    os.system("ls " + data + "/*dab | perl -pe 's/.*\/(.*)\.q?dab/$.\t$+/' > "
              + working + "/datasets.txt")
    # Make networks.bin file
    cmdline = counter + " -k " + working + '/counts/ -o ' + working + \
        '/networks.bin -s ' + working + '/datasets.txt -b ' + working + \
        '/global.txt'
    # Regularization
    if alphas is not None:
        cmdline = cmdline + " -a " + alphas + " -p " + str(pseudo)
    if stddir is not None:
        name = '$+'
        if weights_suffix is not None:
            name += weights_suffix
        # assuming either all dab or all dat
        if os.listdir(stddir)[0].endswith('.dab'):
            os.system("ls " + stddir +
                      "/*\.dab | perl -pe 's/.*\/(.*)\.dab/$.\t" +
                      name + "\t" + name + "/' > " + working + "/contexts.txt")
        else:
            os.system("ls " + stddir +
                      "/*\.dat | perl -pe 's/.*\/(.*)\.dat/$.\t" +
                      name + "\t" + name + "/' > " + working + "/contexts.txt")
    cmdline = cmdline + " -X " + working + "/contexts.txt"
    job.set_name_command(job_name + '-NetBin', cmdline)
    networks_job = job.submit(working + '/NetBin.sh')
    job.set_depends(None)
    return [networks_job, ]

# Counter Predict


def predict(job=None, job_name=None, counter=None, data=None, working=None,
            genome=None, zeros=None, standards=None, threads=16, depends=None,
            weights_suffix=None, pred_queue=None):
    if depends is not None:
        job.set_depends(depends[:])
    predict_jobs = []
    # Make global predictions
    try:
        os.mkdir(working + '/predictions')
    except OSError:
        pass

    cmdline = counter + " -n " + working + '/networks.bin -s ' + working + \
        '/datasets.txt -d ' + data + ' -e ' + genome + ' -o ' + \
        working + '/predictions' + ' -Z ' + zeros

    # run context predict
    if standards is not None:
        contexts_submitted = 0
        # ctxt_jobs = []
        while contexts_submitted < len(standards):
            job_ctxts = \
                standards[contexts_submitted:(contexts_submitted + threads)]
            if weights_suffix is not None:
                job_ctxts = [x + weights_suffix for x in job_ctxts]
            job_thds = len(job_ctxts)
            job_cmd = cmdline + ' -t ' + str(job_thds) + ' ' + \
                ' '.join(job_ctxts)
            job.ppn = job_thds
            job.set_name_command(job_name + '-CtxtPred', job_cmd)
            if pred_queue is not None:
                job.set_queue(pred_queue)
            predict_jobs.append(job.submit(os.path.join(working, 'CtxtPred' +
                                                        str(contexts_submitted)
                                                        + '.sh')))
            contexts_submitted += job_thds
            job.ppn = 1
    job.set_depends(None)
    return predict_jobs

# DChecker Wrapper


def dcheck(job=None, job_name=None, holdout=None, dchecker=None, answers=None,
           working=None, standards=None, stddir=None, depends=None):
    if depends is not None:
        job.set_depends(depends[:])
    dcheck_jobs = []
    try:
        os.mkdir(working + '/dcheck')
    except OSError:
        pass
    # This change is because the DISCOVERY cluster at Dartmouth has trouble
    # with large numbers of DCheck jobs. We can also make this a command array
    # with a minor rewrite if it works better for another cluster.
    dchecks_str = ""
    for context in standards:
        job_cmd = dchecker + ' -w ' + os.path.join(stddir, context + '.dab') \
            + ' -i ' + working + '/predictions/' + context + '.dab '
        job_cmd = dchecker + ' -w ' + os.path.join(stddir, context + '.dab') \
            + ' -i ' + working + '/predictions/' + context + '.dab '
        if holdout is not None:
            job_cmd += ' -g ' + holdout
        job_cmd += ' > ' + working + '/dcheck/' + context + '.auc'
        dchecks_str += job_cmd + '\n'
    # nasty hack for DISCOVERY, writes all commands to one pbs script
    job.set_name_command(job_name + '-DChk', dchecks_str)
    dcheck_jobs.append(job.submit(os.path.join(working, 'Dchk.sh')))
    job.set_depends(None)
    return dcheck_jobs


# CMD LINE PROCESSING
usage = "usage: %prog [options]"
parser = OptionParser(usage, version="%prog dev-0.0.1")
# CONTROL PARAMETERS
# Cross validation
parser.add_option("-S", "--split-gold-standard",
                        dest="split",
                        action="store_true",
                        default=False,
                        help="split gold standard? if not passed, assumed gold \
                        standard is already split by previous run or k=1")
parser.add_option("-L", "--learn-stage",
                        dest="learn",
                        action="store_true",
                        default=False,
                        help="Should the learning stage be run? Depends on \
                        split if cross validation intervals > 1.")
parser.add_option("-N", "--network-stage",
                        dest="network",
                        action="store_true",
                        default=False,
                        help="Should the networks stage be run? Networks \
                        depends on learing having been run.")
parser.add_option("-P", "--predict-stage",
                        dest="predict",
                        action="store_true",
                        default=False,
                        help="Should predict be run? Predict depends on \
                        Networks having been run.")
parser.add_option("-K", "--dcheck",
                        dest="dcheck",
                        help="Should contexts be DChecked after integration? \
                        DCheck depends on Predict having been run.",
                        action="store_true",
                        default=False)

# Queue details
parser.add_option("-Q", "--queue",
                        dest="queue",
                        type="string",
                        default="discovery",
                        help="What cluster parameters should be used \
                        (e.g. discovery, cetus)?")
parser.add_option("-q", "--pred_queue",
                        dest="pred_queue",
                        type="string",
                        help="Which cluster queue to submit prediction jobs \
                        to (adjustable for ctxtpred that takes longer \
                        than 1 day).")
parser.add_option("-E", "--email",
                        dest="email",
                        type="string",
                        help="An e-mail address is required for submission to \
                         some clusters (e.g. discovery).")
# Where do programs exist
parser.add_option("-B", "--sleipnir-binaries-dir",
                        dest="sleipnir",
                        help="DIR containing sleipnir binaries. If not passed, \
                        binaries used must be available in the path.",
                        metavar="DIR")

# INTEGRATION OPTIONS
# Choose Carefully
parser.add_option("-a", "--answers-file",
                        dest="answers",
                        help="FILE with global gold standard relationships -- \
                        only used for gene list.",
                        metavar="FILE")
parser.add_option("-z", "--zeros-file",
                        dest="zeros",
                        help="Zeros file",
                        metavar="FILE")
parser.add_option("-d", "--data-directory",
                        dest="data",
                        help="Directory where data .qdabs are located.",
                        metavar="string")
parser.add_option("-e", "--genes-file",
                        dest="genome",
                        help="Tab-delimited text file containing two "
                        "columns, the first a one-based integer index and "
                        "the second the unique identifier of each gene "
                        "in the genome.",
                        metavar="string")
parser.add_option("-s", "--standards-directory",
                        dest="stddir",
                        help="Directory where standards files are located",
                        metavar="DIR")
parser.add_option("-r", "--alphas-file",
                        dest="alphas",
                        help="Alphas file for bayesian regularization",
                        metavar="FILE")
parser.add_option("-p", "--pseudo-count",
                        dest="pseudo",
                        help="Pseudo count to be used with alphas file "
                        "(default 1)",
                        type="int",
                        default=1,
                        metavar="int")
parser.add_option("-W", "--weights-directory",
                        dest="wtsdir",
                        help="Directory where weights files are located "
                        "(if included, will run weighted integration)",
                        metavar="FILE")
parser.add_option("--weights-suffix",
                  dest="weights_suffix",
                  help="Suffix on weight files (default='_weights') "
                  "note that weight files must have file prefixes "
                  "that match standards files (w/o .dat), o/w can't find them, \
                  avoid periods - confuses Counter",
                  default="_weights",
                  metavar="str")
parser.add_option("-F", "--flipneg-flag",
                        dest="flipneg",
                        help="Flip weights (one minus original) for negative \
                        standards (default=on)",
                        action="store_false", default=True)
parser.add_option("-n", "--noweightneg-flag",
                        dest="noweightneg",
                        help="Use weight one for all negative standards (when \
                        using weight file, default=off)",
                        action="store_true", default=False)
parser.add_option("-m", "--multiplier",
                        dest="multiplier",
                        help="multiplier used for weighting (default = 1)",
                        type="int", default=1)
# Minor Details
parser.add_option("-w", "--working-directory",
                        dest="workdir",
                        help="Perform integration in DIR.",
                        type="string",
                        metavar="DIR")
parser.add_option("-k", "--k-intervals",
                        dest="intervals",
                        type="int",
                        default=1,
                        help="Number of cross validation intervals? 1 = no \
                        cross validation.")
(options, args) = parser.parse_args()

# You always need a gold standard, and CV has to be 1 or greater
if options.answers is None:
    sys.stderr.write("--answers-file file is required.\n")
    sys.exit()
if options.workdir is None:
    sys.stderr.write("--working-directory is required.\n")
    sys.exit()
if options.intervals <= 0:
    sys.stderr.write("--k-intervals must be 1 or greater.\n")
    sys.exit()
if '.' in options.weights_suffix[:-4]:
    sys.stderr.write("--weights_suffix has a period; warning: \
            may mess up networks.bin generation, counter replaces \
            all '.' with '_'\n")


# Add paths if required.
if options.sleipnir is not None:
    for tool in binaries:
        binaries[tool] = options.sleipnir + '/' + binaries[tool]

intervals = []
int_dir = os.path.join(options.workdir, 'intervals')
# Make CV gene sets
if options.intervals > 1:
    intervals = list(range(options.intervals))
    if not os.path.exists(int_dir):
        if not options.split:
            sys.stderr.write('-S must be passed if intervals have not already \
            been created.\n')
            sys.exit()
        else:
            os.mkdir(int_dir)
    if options.split:
        os.system(binaries['Dat2Dab'] + ' -E -i ' + options.answers + ' > ' +
                  os.path.join(int_dir, 'all.txt'))
        all_genes = set()
        all_file = open(os.path.join(int_dir, 'all.txt'))
        for line in all_file:
            all_genes.add(line.strip())
        all_file.close()
        all_genes = list(all_genes)
        per_each = len(all_genes) // options.intervals
        extra_genes = len(all_genes) % options.intervals
        random.shuffle(all_genes)
        used_already = 0
        for interval in intervals:
            in_interval = per_each
            if interval < extra_genes:
                in_interval += 1
            interval_genes = \
                all_genes[used_already:(used_already + in_interval)]
            used_already += in_interval
            interval_file = open(os.path.join(int_dir,
                                              str(interval) + '.txt'), 'w')
            interval_file.write('\n'.join(interval_genes) + '\n')
            interval_file.close()

# Make directories for each interval
if not os.path.exists(os.path.join(options.workdir, 'all')):
    os.mkdir(os.path.join(options.workdir, 'all'))
for interval in intervals:
    if not os.path.exists(os.path.join(options.workdir, str(interval))):
        os.mkdir(os.path.join(options.workdir, str(interval)))

# Find standards
standards = []
potential = os.listdir(options.stddir)
for item in potential:
    if item.endswith('.dab') or item.endswith('.dat'):
        standards.append(item.split('.')[0])

# Find weight files if they exist
weights = {}
if options.wtsdir:
    for std in standards:
        wfile_path = os.path.join(options.wtsdir, std + options.weights_suffix)
        if os.path.exists(wfile_path + '.dat'):
            wfile_path += '.dat'
        elif os.path.exists(wfile_path + '.dab'):
            wfile_path += '.dab'
        else:
            sys.stderr.write('Could not find corresponding weights file for '
                             + std + '\n')
            sys.exit()

        weights[std] = wfile_path
else:
    options.weights_suffix = None


job, local_cmds = None, None
if options.queue == 'discovery':
    from pbsjob import PBSJob
    job = PBSJob(addr=options.email,
                 command="echo test",
                 walltime="47:59:00",
                 queue="largeq")
elif options.queue == 'slurm':
    from slurmjob import SLURMJob
    job = SLURMJob(addr=options.email,
                   command="echo test",
                   walltime="24:00:00",
                   queue="1day",
                   workdir=options.workdir)
else:
    local_cmds = defaultdict(str)

# Check requirements
if options.learn or options.network or options.predict:
    if options.data is None:
        sys.stderr.write("--data-directory is required.\n")
        sys.exit()
if options.predict:
    if options.genome is None:
        sys.stderr.write("--genes-file is required.\n")
        sys.exit()

# Run non-cv
depend_jobs = None
if options.learn:
    depend_jobs = learn(job=job, job_name='all', counter=binaries['Counter'],
                        global_answers=options.answers, data=options.data,
                        working=os.path.join(options.workdir, 'all'),
                        zeros=options.zeros, standards=standards,
                        stddir=options.stddir, wtsdir=options.wtsdir,
                        weights=weights, flipneg=options.flipneg,
                        multiplier=options.multiplier,
                        noweightneg=options.noweightneg,
                        local_cmds=local_cmds)
if options.network:
    depend_jobs = networks(job=job, job_name='all',
                           counter=binaries['Counter'], data=options.data,
                           working=os.path.join(options.workdir, 'all'),
                           alphas=options.alphas, pseudo=options.pseudo,
                           stddir=options.stddir, depends=depend_jobs,
                           weights_suffix=options.weights_suffix)
if options.predict:
    depend_jobs = predict(job=job, job_name='all', counter=binaries['Counter'],
                          data=options.data,
                          working=os.path.join(options.workdir, 'all'),
                          genome=options.genome, zeros=options.zeros,
                          standards=standards, depends=depend_jobs,
                          weights_suffix=options.weights_suffix,
                          pred_queue=options.pred_queue)
if options.dcheck:
    depend_jobs = dcheck(job=job, job_name='all',
                         dchecker=binaries['DChecker'],
                         answers=options.answers,
                         working=os.path.join(options.workdir, 'all'),
                         standards=standards, stddir=options.stddir,
                         depends=depend_jobs)


print(list(local_cmds.values())[0])
os.system(list(local_cmds.values())[0])


for interval in intervals:
    job_name = 'cv' + str(interval)
    holdout = os.path.join(int_dir, str(interval) + '.txt')
    depend_jobs = None
    if options.learn:
        depend_jobs = learn(job=job, job_name=job_name, holdout=holdout,
                            counter=binaries['Counter'],
                            global_answers=options.answers,
                            data=options.data,
                            working=os.path.join(options.workdir,
                                                 str(interval)),
                            zeros=options.zeros, standards=standards,
                            stddir=options.stddir)
    if options.network:
        depend_jobs = networks(job=job, job_name=job_name,
                               counter=binaries['Counter'], data=options.data,
                               working=os.path.join(options.workdir,
                                                    str(interval)),
                               alphas=options.alphas, pseudo=options.pseudo,
                               stddir=options.stddir, depends=depend_jobs)
    if options.predict:
        depend_jobs = predict(job=job, job_name=job_name,
                              counter=binaries['Counter'], data=options.data,
                              working=os.path.join(options.workdir,
                                                   str(interval)),
                              genome=options.genome, zeros=options.zeros,
                              standards=standards, depends=depend_jobs)
    if options.dcheck:
        depend_jobs = dcheck(job=job, job_name=job_name, holdout=holdout,
                             dchecker=binaries['DChecker'],
                             answers=options.answers,
                             working=os.path.join(options.workdir,
                                                  str(interval)),
                             standards=standards, stddir=options.stddir,
                             depends=depend_jobs)
