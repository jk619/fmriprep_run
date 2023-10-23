"""
Heuristics file in use at the NYU's Center for Brain Imaging
(Our data come from a Siemens Prisma 3T scanner).
At the end of the classification, the function
`clean_up_unneeded_tags_from_info` simplifies the BIDS file
names

Author: Pablo Velasco
Date: 4/28/2020
"""

import os

# Note: When we have pairs of scans "original" + "normalize", we only
#       keep the normalized image in the BIDS structure. The DICOMs for
#       the original images are extracted in the 'sourcedata' folder.

DEFAULT_FIELDS = {
    # For the `dataset_description.json`:
    "Acknowledgements":
        "We thank Pablo Velasco and the rest of the NYU CBI (Center "
        "for Brain Imaging) personnel for preparing the BIDS "
        "dataset. "
        "TODO: whom else you want to acknowledge"
}


def create_key(subdir, file_suffix, outtype=('nii.gz', 'dicom'),
               annotation_classes=None, prefix=''):
    if not subdir:
        raise ValueError('subdir must be a valid format string')
    # may be even add "performing physician" if defined??
    template = os.path.join(
        prefix,
        "{bids_subject_session_dir}",
        subdir,
        "{bids_subject_session_prefix}_%s" % file_suffix
    )
    return template, outtype, annotation_classes


def find_PE_direction_from_protocol_name(prot_name, default_dir_name='normal'):
    # valid phase-encoding directions in the protocol name
    PE_directions = ['AP','PA','RL','LR','rev']
    direction = default_dir_name
    for peDir in PE_directions:
        if (
            '_'+peDir in prot_name
            or '-'+peDir in prot_name
        ):
            direction = peDir
            break

    return direction


def extract_task_name(prot_name):
    """ extracts task name from the protocol_name """

    prot_name = prot_name.lower()
    known_tasks=[
        'rest',
        'face',
        'conflict',
        'gamble',
        'fix',
        'films',
        'inscapes',
    ]
    for task_name in known_tasks:
        if task_name in prot_name:
            break    # we keep the current value for "task"
        else:
            task_name = ''

    # if we don't find a known task, try finding the keyword "task":
    if task_name == '':
        if 'task' in prot_name and not prot_name.endswith('task'):
            # we want to capture what comes after "task", up to the next
            #    dash ("-") or underscore ("_"):
            task_name = prot_name.split('task')[1]
            # remove any initial "-" or "_":
            if task_name[0] == '-' or task_name[0] == '_':
                task_name = task_name[1:]
            # discard anything after the next "-" or "_":
            task_name = task_name.split('_')[0]
            task_name = task_name.split('-')[0]
            # remove any spaces we might have:
            task_name = task_name.replace(" ", "")
        else:
            task_name = 'TASK'    # fallback.  BIDS requires an alphanumeric string (no spaces...)

    return task_name


def add_series_to_info_dict(series_id, mykey, info, acq=''):
    """ adds a series to the 'info' dictionary """

    if info == None or mykey == '':
        return error

    if mykey in info:
        if acq == '':
            info[mykey].append({'item': series_id})
        else:
            info[mykey].append({'item': series_id, 'acq': acq})
    else:
        # if it isn't, add this key, specifying the first item:
        if acq == '':
            info[mykey] = [{'item': series_id}]
        else:
            info[mykey] = [{'item': series_id, 'acq': acq}]


def previous_series_is_sbref(seqinfo, index):
    """
    Checks if the previous series is a SBRef image
    seqinfo: list with the sequences information
    index: index of the series, so the previous series would be index-1
    """
    return (
        index > 0
        and seqinfo[index - 1].series_description.lower().endswith('_sbref')
    )


def another_series_has_identical_start_time(seqinfo, index):
    """
    Check if check if another series has identical start date and
    time.
    If so, they are just the same acquisition and the only difference
    is that one of them would be intensity-normalized or motion-
    corrected, etc.
    """

    if not hasattr(seqinfo[index],'time'):
        return False
    else:
        series_with_identical_start_time = [
            s for s in seqinfo if (
                s.date and s.date==seqinfo[index].date
                and s.time and s.time==seqinfo[index].time
            )
        ]
        return len(series_with_identical_start_time)>1


def clean_up_unneeded_tags_from_info(info):
    """
    For "info" keys that only have one run, delete the "_run-<run>" tag
    """
    # Because we are going to be modifying the keys of the dictionary,
    # we get the original values now:
    orig_keys = [k for k in info.keys()]
    for k in orig_keys:
        if (
            len(info[k]) == 1
            and '_run-{item:02d}' in k[0]
        ):
            # For entries in which run number is not hard-coded:
            new_k = list(k)
            new_k[0] = new_k[0].replace('_run-{item:02d}','')
            # remove the old k from "info" and replace it with new_k
            # while keeping the same value
            info[tuple(new_k)] = info.pop(k)
        elif '_run-01' in k[0]:
            # If "_run-01" is hardcoded, check to see if there is a
            # corresponding "_run-02". If not, remove "_run-01" from the
            # key:
            run2 = k[0].replace('_run-01','_run-02')
            if not run2 in [kk[0] for kk in info.keys()]:
                new_k = list(k)
                new_k[0] = new_k[0].replace('_run-01','')
                # remove the old k from "info" and replace it with new_k
                # while keeping the same value
                info[tuple(new_k)] = info.pop(k)

    # Now, go through the "func" runs and if all of them have the same
    # acquisition, remove the "_acq-<acq>" tag from the keys:
    func_keys = [k for k in info.keys() if '/func/' in k[0]]
    # use a "set" to keep only unique values:
    acqs = set([info[k][0]['acq'] for k in func_keys])
    if len(acqs) == 1:
        for k in func_keys:
            new_k = list(k)
            new_k[0] = new_k[0].replace('_acq-{acq}','')
            info[tuple(new_k)] = info.pop(k)


def sort_seqinfo_by_series_date_time(seqinfo):
    """
    Sorts the 'seqinfo' list according to the series acquisition date
    and time (using the "series_uid" to make sure there are no repeats)

    This assures that all the entries are in acquisition order (which
    is not always the case).

    input: seqinfo: original seqinfo list
    output: sortedSeqinfo: seqinfo sorted list
    """
    # sort by concatenated date and time (strings):
    # (use set to only keep unique ones):
    dateTimes = set([s.date+s.time for s in seqinfo if s.date and s.time])
    sortedSeqinfo = []
    for dt in sorted(dateTimes):
        dTseries = [
            ss for ss in seqinfo if (
                ss.date
                and ss.time
                and ss.date+ss.time == dt
            )
        ]
        # sort series with identical date and time by series_uid:
        for sid in sorted([s.series_uid for s in dTseries if s.series_uid]):
            for ss in dTseries:
                if ss.series_uid == sid:
                    sortedSeqinfo.append(ss)
                    # add the following series if they don't have date and time:
                    idx = seqinfo.index(ss) + 1
                    while (
                            idx < len(seqinfo)
                            and not (seqinfo[idx].date and seqinfo[idx].time)
                    ):
                        sortedSeqinfo.append(seqinfo[idx])
                        idx += 1

        # Now, add the series which do not have series_uid:
        for ss in dTseries:
            if ss not in sortedSeqinfo:
                sortedSeqinfo.append(ss)

    # Now, add any missing series:
    for ss in seqinfo:
        if ss not in sortedSeqinfo:
            sortedSeqinfo.append(ss)

    return sortedSeqinfo


def infotodict(seqinfo):
    """
    Heuristic evaluator for determining which runs belong where allowed
    template fields - follow python string module:

    item: index within category
    subject: participant id
    seqitem: run number during scanning
    subindex: sub index within group
    """
    ###   NIFTI and DICOM   ###
    # anat:
    t1 = create_key('anat','acq-{acq}_run-{item:02d}_T1w')
    t2 = create_key('anat','acq-{acq}_run-{item:02d}_T2w')
    pd_bias_body = create_key('anat','acq-biasBody_run-{item:02d}_PD')
    pd_bias_receive = create_key('anat','acq-biasReceive_run-{item:02d}_PD')

    # func:
    # For the functional runs, we want to have a different key for each
    # task. Since we don't know a-priori what the task names will be, we
    # don't create the different keys here, but we'll create them
    # dynamically, as needed.

    # diffusion:
    dwi = create_key('dwi','acq-{acq}_run-{item:02d}_dwi')
    dwi_sbref = create_key('dwi','acq-{acq}_run-{item:02d}_sbref')

    # fmaps:
    # For the SE-EPI field-map runs (for topup), both for fmri and dwi,
    # we want to have a different key for each PE direction.  Rather
    # than creating the different keys here, we'll create them
    # dynamically, as needed.
    fmap_gre_mag = create_key('fmap','acq-GRE_run-{item:02d}_magnitude')       # GRE fmap
    fmap_gre_phase = create_key('fmap','acq-GRE_run-{item:02d}_phasediff')     #

    ###  DICOM only   ###
    # These are images we don't want for analysis, but we still want to
    # keep a copy of the DICOMs in the 'sourcedata' folder.  We manage
    # that by specifying the 'outtype' to be only 'dicom':
    # anat:
    t1_scout = create_key('anat','acq-Scout_run-{item:02d}_T1w', outtype = ('dicom',))
    t1_dicom = create_key('anat','acq-{acq}_run-{item:02d}_T1w', outtype = ('dicom',))
    t2_dicom = create_key('anat','acq-{acq}_run-{item:02d}_T2w', outtype = ('dicom',))
    # Misc:
    phoenix_doc = create_key('misc','phoenix', outtype = ('dicom',))

    info = {t1:[], t2:[], pd_bias_body:[], pd_bias_receive:[],
            dwi:[], dwi_sbref:[],
            fmap_gre_mag:[], fmap_gre_phase:[],
            t1_scout:[], t1_dicom: [], t2_dicom: [], phoenix_doc:[]}
    run_numbers = {}    # dict with the run number for each task

    # because not always are series sorted by acquisition date/time,
    # sort them:
    seqinfo = sort_seqinfo_by_series_date_time(seqinfo)

    for idx, s in enumerate(seqinfo):
        # s is a namedtuple with fields equal to the names of the columns
        # found in the dicominfo.tsv file

        # check the PE direction:
        direction = find_PE_direction_from_protocol_name(
            s.protocol_name, default_dir_name=''
        )
        # we don't need the original case of the protocol name:
        prot_name = s.protocol_name.lower()
        acq = ''

        ###   T1w   ###
        # 1) Auto-Align scout (the original images, not the derived):
        #    look for "AA*Scout" or "AA*scout" in protocol_name:
        if (
            'scout' in prot_name
            and s.is_derived == False
        ):
            info[t1_scout].append({'item': s.series_id})

        # 2) High resolution T1w:
        # single volume, protocol name including T1, MPRAGE, MP-RAGE, MPR,...
        elif (
            s.dim4 == 1
            and 'fl' in s.sequence_name
            and (
                't1' in prot_name
                or ('mp' in prot_name and 'rage' in prot_name)
                or 'mpr' in prot_name
            )
        ):
            acq = 'highres' + direction   # if direction is empty, aqc='highres'

            # If this image is NOT normalized, check if another
            # series has identical acquisition date and time.
            # If so, we'll keep only the normalized version, and for
            # this one we keep just the DICOM.
            # (Note: older versions of heudiconv don't include 'time'):
            # Otherwise, we extract it:
            if (
                'NORM' not in s.image_type
                and another_series_has_identical_start_time(seqinfo, idx)
            ):
                info[t1_dicom].append({'item': s.series_id, 'acq': acq})
            else:
                info[t1].append({'item': s.series_id, 'acq': acq})

        # 3) FSE T1w:
        # single volume, series description includes TSE or FSE,
        # protocol name includes T1
        elif (
            s.dim4 == 1
            and 't1' in prot_name
            and 'tse' in s.sequence_name
        ):
            acq = 'fse' + direction   # if direction is empty, aqc='fse'
            info[t1].append({'item': s.series_id, 'acq': acq})

        # 4) ME-MPRAGE:
        # several volumes, series description includes "tfl_me3d":
        elif (
            s.sequence_name.startswith('tfl_me3d')
        ):
            acq = 'highres' + direction   # if direction is empty, aqc='highres'
            info[t1].append({'item': s.series_id, 'acq': acq})


        ###   T2w   ###
        # 1) Standard high-res T2w used for cortical segmentation:
        # single volume, protocol name including T2, T2w, TSE, SPACE, SPC:
        elif (
            s.dim4 == 1
            and (
                'T2' in s.protocol_name
                or 'tse' in prot_name
                or 'space' in prot_name
                or 'spc' in prot_name
            )
        ):
            acq = 'highres' + direction   # note: if direction is empty, aqc='highres'

            # If this image is NOT normalized, check if another
            # series has identical acquisition date and time.
            # If so, we'll keep only the normalized version, and for
            # this one we keep just the DICOM.
            # (Note: older versions of heudiconv don't include 'time'):
            # Otherwise, we extract it:
            if (
                'NORM' not in s.image_type
                and another_series_has_identical_start_time(seqinfo, idx)
            ):
                info[t2_dicom].append({'item': s.series_id, 'acq': acq})
            else:
                info[t2].append({'item': s.series_id, 'acq': acq})

        # 2) Fast Spin-Echo used as a localizer:
        # single volume, sequence name: 'h2d1' ('haste')"
        elif (
            s.dim4 == 1
            and 'h2d1' in s.sequence_name
        ):
            info[t2_dicom].append({'item': s.series_id, 'acq': 'haste'})

        ###   PD   ###
        # BIAS images  (for coil sensitivity estimation) are typically
        # PD-weighted
        elif (
            'bias' in prot_name
            and 'tfl3d' in s.sequence_name
        ):
            if 'body' in prot_name:
                info[pd_bias_body].append({'item': s.series_id})
            else:
                info[pd_bias_receive].append({'item': s.series_id})

        ###   FUNCTIONAL   ###
        # We want to make sure the _SBRef, PhysioLog and phase series
        # (if present) are labeled the same as the main (magnitude)
        # image. So we only focus on the magnitude series (to exclude
        # phase images) and more than 3 volumes (to exclude _SBRef) and
        # then we search if the phase and/or _SBRef are present.
        elif (
            s.dim4 >= 4
            and 'epfid2d' in s.sequence_name
            and (
                'M' in s.image_type
                or 'FMRI' in s.image_type
            )
            and not s.series_description.lower().endswith('_sbref')
            and not 'DERIVED' in s.image_type
        ):

            # check PE direction:
            # we will write the direction under the 'acq-' tag.
            acq = direction or 'normal'

            # check task name:
            task = extract_task_name(s.protocol_name)
            # find the run number for this task:
            if run_numbers.get(task):
                run_numbers[task] += 1
            else:
                run_numbers[task] = 1

            # dictionary key specific for this task type:
            mykey = create_key(
                'func',
                'task-%s_acq-{acq}_run-%02d_bold' % (task, run_numbers[task])
            )
            add_series_to_info_dict(s.series_id, mykey, info, acq)
            next_series = idx+1    # used for physio log below

            ###   is phase image present?   ###
            # At least for Siemens systems, if magnitude/phase was
            # selected, the phase images come as a separate series
            # immediatelly following the magnitude series.
            # (note: make sure you don't check beyond the number of
            # elements in seqinfo...)
            if (
                idx+1 < len(seqinfo)
                and 'P' in seqinfo[idx+1].image_type
            ):
                mykey_pha = create_key(
                    'func',
                    'task-%s_acq-{acq}_run-%02d_phase' % (task, run_numbers[task])
                )
                add_series_to_info_dict(
                    seqinfo[idx + 1].series_id, mykey_pha, info, acq
                )
                next_series = idx+2    # used for physio log below

            ###   SB REF   ###
            # here, within the functional run code, check to see if the
            #  previous run protocol name ended in _SBREF, to assign the
            #  same task name and run number.
            if previous_series_is_sbref(seqinfo, idx):
                mykey_sb = create_key(
                    'func',
                    'task-%s_acq-{acq}_run-%02d_sbref' % (task, run_numbers[task])
                )
                add_series_to_info_dict(
                    seqinfo[idx - 1].series_id, mykey_sb, info, acq
                )

            ###   PHYSIO LOG   ###
            # here, within the functional run code, check to see if the
            #  next run image_type lists "PHYSIO", to assign the
            #  same task name and run number.
            if (
                next_series < len(seqinfo)
                and 'PHYSIO' in seqinfo[next_series].image_type
            ):
                mykey_physio = create_key(
                    'func',
                    'task-%s_acq-{acq}_run-%02d_physio' % (task, run_numbers[task]),
                    outtype = ('dicom','physio')
                )
                add_series_to_info_dict(
                    seqinfo[next_series].series_id, mykey_physio, info, acq
                )


        ###   FIELD MAPS   ###

        # A) Spin-Echo distortion maps to be used with topup:
        #  (note: we only look for magnitude images, because we'll catch
        #   the phase images below
        #   Also, we exclude _sbref because the sbref from MB diffusion
        #   also have epse2d in the sequence_name
        elif (
            s.dim4 <= 3
            and 'epse2d' in s.sequence_name
            and 'M' in s.image_type
            and (
                'dist' in prot_name
                or 'map' in prot_name
                or 'field' in prot_name
            )
            and not s.series_description.lower().endswith('_sbref')
        ):
            # check PE direction:
            run_dir = direction or 'normal'

            # dictionary key specific for this SE-fmap direction:
            mykey = create_key(
                'fmap',
                'acq-fMRI_dir-%s_run-{item:02d}_epi' % run_dir
            )
            add_series_to_info_dict(s.series_id, mykey, info)

        # B) GRE fmap:
        elif (
            'fm2d' in s.sequence_name
            and ('fieldmap' in prot_name or 'field_map' in prot_name)
        ):
            if s.image_type[2] == 'M':
                # magnitude image
                info[fmap_gre_mag].append({'item': s.series_id})
            elif s.image_type[2] == 'P':
                # phase image:
                info[fmap_gre_phase].append({'item': s.series_id})


        ###   DIFFUSION   ###

        # We could also check: (s.image_type[2] == 'DIFFUSION')
        elif 'ep_b' in s.sequence_name:
            # Siemens product diffusion sequence

            if 'fwf' in prot_name.lower():
                for acq in ('lte', 'pte', 'ste'):
                    if '_' + acq in prot_name.lower():
                        if '_PA' in prot_name:
                            acq = acq + '_dir-PA'
                        # dictionary key specific for this acquisition type:
                        mykey = create_key(
                            'dwi',
                            'acq-%s_run-{item:02d}_epi' % acq
                        )
                        add_series_to_info_dict(s.series_id, mykey, info)

            # This is not very rigorous, but in general, diffusion runs
            # will have more than a couple of volumes.  I'll use 5, to
            # be safe:
            elif s.dim4 >= 5:
                # this is a standard diffusion acquisition
                acq = str(s.dim4)+'vols'
                info[dwi].append({'item': s.series_id, 'acq': acq})

                # check to see if the previous run is a SBREF:
                if (
                    previous_series_is_sbref(seqinfo, idx)
                    and 'epse2d' in seqinfo[idx - 1].sequence_name
                ):
                    info[dwi_sbref].append(
                        {'item': seqinfo[idx - 1].series_id, 'acq': acq}
                    )

            else:
                # this is a fmap for diffusion.

                # Get the PE direction, for topup:
                run_dir = direction or 'normal'

                # dictionary key specific for this SE-fmap direction:
                mykey = create_key(
                    'fmap',
                    'acq-dwi_dir-%s_run-{item:02d}_epi' % run_dir
                )
                add_series_to_info_dict(s.series_id, mykey, info)

                if previous_series_is_sbref(seqinfo, idx):
                    # TO-DO: for now, extract the _sbref dwi fmap image
                    # as DICOMs, because BIDS doesn't allow them.
                    mykey = create_key(
                        'fmap',
                        'acq-dwi_dir-%s_run-{item:02d}_sbref' % run_dir,
                        outtype = ('dicom',)
                    )
                    add_series_to_info_dict(
                        seqinfo[idx - 1].series_id, mykey, info
                    )

        ###   PHOENIX FILE   ###

        elif (
            'PhoenixZIPReport' in s.series_description
            and s.image_type[3] == 'CSA REPORT'
        ):
            info[phoenix_doc].append({'item': s.series_id})

    clean_up_unneeded_tags_from_info(info)
    return info
