#!/bin/bash -l
#!/home/amax/anaconda3/envs/mrt3/bin/python3.7 "$@"

import glob, json, math, os, re, shutil, subprocess
from distutils.spawn import find_executable
from mrtrix3 import app, file, fsl, image, path, run


__version__ = 'BIDS-App \'MRtrix3_connectome\' standalone'
option_prefix = '-'

def runSubject(bids_dir, label, output_prefix):

  output_dir = os.path.join(output_prefix, label)
  if os.path.exists(output_dir):
    app.warn('Output directory for subject \'' + label + '\' already exists; contents will be erased when this execution completes')

  fsl_path = os.environ.get('FSLDIR', '')
  if not fsl_path:
    app.error('Environment variable FSLDIR is not set; please run appropriate FSL configuration script')

  robex_found = find_executable('ROBEX') and find_executable('runROBEX.sh')
  N4_found = find_executable('N4BiasFieldCorrection')

  if robex_found and N4_found:
    app.console('N4BiasFieldCorrection and ROBEX found; will use for bias field correction and brain extraction')
    brain_extraction_cmd = 'runROBEX.sh'
  else:
    if robex_found and not N4_found:
      app.console('N4BiasFieldCorrection not found; will use FSL for bias field correction & brain extraction')
    elif N4_found and not robex_found:
      app.console('ROBEX not found; will use FSL for brain extraction')
    else:
      app.console('N4BiasFieldCorrection and ROBEX not found; will use FSL for brain extraction')
    brain_extraction_cmd = fsl.exeName('fsl_anat')

  dwibiascorrect_algo = 'ants' # modified by chenfei
  if not N4_found:
    # Can't use fsl.exeName() here, since we want to proceed even if it's not found
    if find_executable('fast') or find_executable('fsl5.0-fast'):
      dwibiascorrect_algo = 'fsl' # modified by chenfei
      app.console('Could not find ANTs program N4BiasFieldCorrection; '
                  'using FSL FAST for bias field correction')
    else:
      dwibiascorrect_algo = ''
      app.warn('Could not find ANTs program \'N4BiasFieldCorrection\' or FSL program \'fast\'; '
               'will proceed without performing initial DWI bias field correction')

  if not app.args.parcellation:
    app.error('For participant-level analysis, desired parcellation must be provided using the ' + option_prefix + 'parcellation option')

  if app.args.template_reg:
    if app.args.parcellation in [ 'aal', 'aal2', 'craddock200', 'craddock400', 'perry512' ]:
      template_registration_software = app.args.template_reg
    else:
      app.warn('Template registration not being performed; ' + option_prefix + 'template_reg option ignored')
      template_registration_software = ''
  else:
    template_registration_software = 'ants' if app.args.parcellation in [ 'aal', 'aal2', 'craddock200', 'craddock400', 'perry512' ] else ''
  if template_registration_software == 'ants':
    if not find_executable('ANTS') or not find_executable('WarpImageMultiTransform'):
      app.error('Commands \'ANTS\' and \'WarpImageMultiTransform\' must be present in PATH to use ANTs software for template registration')

  template_image_path = ''
  template_mask_path = ''
  parc_image_path = ''
  parc_lut_file = ''
  mrtrix_lut_dir = '/home/amax/anaconda3/share/mrtrix3/labelconvert' #modified by chenfei
  mrtrix_lut_file = ''
  # MNI152NLin2009cAsym_path = os.path.join(bids_dir, 'derivatives', 'atlas', 'tpl-MNI152NLin2009cAsym_res-02_T1w.nii.gz')

  if app.args.parcellation in [ 'desikan', 'destrieux', 'hcpmmp1' ]:
    if not 'FREESURFER_HOME' in os.environ:
      app.error('Environment variable FREESURFER_HOME not set; please verify FreeSurfer installation')
    freesurfer_subjects_dir = os.environ['SUBJECTS_DIR'] if 'SUBJECTS_DIR' in os.environ else os.path.join(os.environ['FREESURFER_HOME'], 'subjects')
    if not os.path.isdir(freesurfer_subjects_dir):
      app.error('Could not find FreeSurfer subjects directory (expected location: ' + freesurfer_subjects_dir + ')')
    for subdir in [ 'fsaverage', 'lh.EC_average', 'rh.EC_average' ]:
      if not os.path.isdir(os.path.join(freesurfer_subjects_dir, subdir)):
        app.error('Could not find requisite FreeSurfer subject directory \'' + subdir + '\' (expected location: ' + os.path.join(freesurfer_subjects_dir, subdir) + ')')
    if app.args.parcellation == 'hcpmmp1':
      if not all([os.path.isfile(os.path.join(freesurfer_subjects_dir, 'fsaverage', 'label', hemi + 'h.HCPMMP1.annot')) for hemi in [ 'l', 'r' ]]):
        app.error('Could not find necessary annotation labels for applying HCPMMP1 parcellation ' + \
                  '(expected location: ' + os.path.join(freesurfer_subjects_dir, 'fsaverage', 'label', '?h.HCPMMP1.annot'))
    reconall_path = find_executable('recon-all')
    if not reconall_path:
      app.error('Could not find FreeSurfer script recon-all; please verify FreeSurfer installation')
    if app.args.parcellation == 'hcpmmp1':
      for cmd in [ 'mri_surf2surf', 'mri_aparc2aseg' ]:
        if not find_executable(cmd):
          app.error('Could not find FreeSurfer command ' + cmd + ' (necessary for applying HCPMMP1 parcellation); please verify FreeSurfer installation')
      if not all(os.path.isfile(os.path.join(freesurfer_subjects_dir, 'fsaverage', 'label', hemi + 'h.HCPMMP1.annot')) for hemi in [ 'l', 'r' ]):
        app.error('Could not find FreeSurfer annotation files for HCPMMP1 parcellation (expected location: ' + os.path.join(freesurfer_subjects_dir, 'fsaverage', 'label', '?h.HCPMMP1.annot'))
    # Query contents of recon-all script, looking for "-openmp" and "-parallel" occurences
    # Add options to end of recon-all -all call, based on which of these options are available
    #   as well as the value of app.numThreads
    # - In 5.3.0, just the -openmp option is available
    # - In 6.0.0, -openmp needs to be preceded by -parallel
    reconall_multithread_options = []
    if app.numThreads is None or app.numThreads:
      with open(reconall_path, 'r') as f:
        reconall_text = f.read().splitlines()
      for line in reconall_text:
        line = line.strip()
        if line == 'case "-parallel":':
          reconall_multithread_options = [ '-parallel' ] + reconall_multithread_options
        # If number of threads in this script is not being explicitly controlled,
        #   allow recon-all to use its own default number of threads
        elif line == 'case "-openmp":' and app.numThreads is not None:
          reconall_multithread_options.extend([ '-openmp', str(app.numThreads) ])
    if reconall_multithread_options:
      reconall_multithread_options = ' ' + ' '.join(reconall_multithread_options)
    else:
      reconall_multithread_options = ''
    app.var(reconall_multithread_options)
    if app.args.parcellation in [ 'desikan', 'destrieux' ]:
      parc_lut_file = os.path.join(os.environ['FREESURFER_HOME'], 'FreeSurferColorLUT.txt')
    else:
      parc_lut_file = os.path.join(mrtrix_lut_dir, 'hcpmmp1_original.txt')
    if app.args.parcellation == 'desikan':
      mrtrix_lut_file = os.path.join(mrtrix_lut_dir, 'fs_default.txt')
    elif app.args.parcellation == 'destrieux':
      mrtrix_lut_file = os.path.join(mrtrix_lut_dir, 'fs_a2009s.txt')
    else:
      mrtrix_lut_file = os.path.join(mrtrix_lut_dir, 'hcpmmp1_ordered.txt')

  elif app.args.parcellation in [ 'aal', 'aal2', 'craddock200', 'craddock400', 'perry512' ]:
    template_image_path = os.path.join(fsl_path, 'data', 'standard', 'MNI152_T1_2mm.nii.gz')
    template_mask_path = os.path.join(fsl_path, 'data', 'standard', 'MNI152_T1_2mm_brain_mask.nii.gz')
    if app.args.parcellation == 'aal':
      parc_image_path = os.path.abspath(os.path.join(os.sep, 'opt', 'aal', 'ROI_MNI_V4.nii'))
      parc_lut_file = os.path.abspath(os.path.join(os.sep, 'opt', 'aal', 'ROI_MNI_V4.txt'))
      mrtrix_lut_file = os.path.join(mrtrix_lut_dir, 'aal.txt')
    elif app.args.parcellation == 'aal2':
      parc_image_path = os.path.abspath(os.path.join(os.sep, 'opt', 'aal', 'ROI_MNI_V5.nii'))
      parc_lut_file = os.path.abspath(os.path.join(os.sep, 'opt', 'aal', 'ROI_MNI_V5.txt'))
      mrtrix_lut_file = os.path.join(mrtrix_lut_dir, 'aal2.txt')
    else:
      parc_lut_file = mrtrix_lut_file = ''
      if app.args.parcellation == 'craddock200':
        parc_image_path = os.path.abspath(os.path.join(os.sep, 'opt', 'ADHD200_parcellate_200.nii.gz'))
      elif app.args.parcellation == 'craddock400':
        parc_image_path = os.path.abspath(os.path.join(os.sep, 'opt', 'ADHD200_parcellate_400.nii.gz'))
      elif app.args.parcellation == 'perry512':
        parc_image_path = os.path.abspath(os.path.join(os.sep, 'opt', '512inMNI.nii'))

  def findAtlasFile(filepath, description):
    if not filepath:
      return ''
    if os.path.isfile(filepath):
      return filepath
    if not getattr(app.args, 'atlas_path', None):
      app.error('Could not find ' + description + ' (expected location: ' + filepath + ')')
    newpath = os.path.join(os.path.dirname(app.args.atlas_path), os.path.basename(filepath))
    if os.path.isfile(newpath):
      return newpath
    app.error('Could not find ' + description + ' (tested locations: \'' + filepath + '\', \'' + newpath + '\')')
    return None

  template_image_path = findAtlasFile(template_image_path, 'template image')
  template_mask_path = findAtlasFile(template_mask_path, 'template brain mask image')
  parc_image_path = findAtlasFile(parc_image_path, 'parcellation image')
  parc_lut_file = findAtlasFile(parc_lut_file, 'parcellation lookup table file')


  if mrtrix_lut_file and not os.path.exists(mrtrix_lut_file):
    app.error('Could not find MRtrix3 connectome lookup table file (expected location: ' + mrtrix_lut_file + ')')

  app.makeTempDir()

  app.console('Importing DWI data into temporary directory')
  dwi_image_list = glob.glob(os.path.join(bids_dir, label, 'dwi', label + '*_dwi.nii*'))
  dwi_index = 1
  for entry in dwi_image_list:
    # os.path.split() falls over with .nii.gz extensions; only removes the .gz
    prefix = entry.split(os.extsep)[0]
    if os.path.isfile(prefix + '.bval') and os.path.isfile(prefix + '.bvec'):
      prefix = prefix + '.'
    else:
      prefix = os.path.join(bids_dir, 'dwi')
      if not (os.path.isfile(prefix + 'bval') and os.path.isfile(prefix + 'bvec')):
        app.error('Unable to locate valid diffusion gradient table for image \'' + entry + '\'')
    grad_import_option = ' -fslgrad ' + prefix + 'bvec ' + prefix + 'bval'
    json_path = prefix + 'json'
    json_import_option = ''
    if os.path.isfile(json_path):
      json_import_option = ' -json_import ' + json_path
    elif not app.args.preprocessed:
      app.error('No sidecar JSON file found for image \'' + entry + '\'; cannot proceed with DWI preprocessing without this information')
    run.command('mrconvert ' + entry + grad_import_option + json_import_option + ' ' + path.toTemp('dwi' + str(dwi_index) + '.mif', True))
    shutil.copy(json_path, path.toTemp('dwi' + str(dwi_index) + '.json', True))
    dwi_index += 1

  # Go hunting for reversed phase-encode data dedicated to field map estimation
  fmap_image_list = []
  fmap_dir = os.path.join(bids_dir, label, 'fmap')
  fmap_index = 1
  if os.path.isdir(fmap_dir):
    if app.args.preprocessed:
      app.error('fmap/ directory detected for subject \'' + label + '\' despite use of ' + option_prefix + 'preprocessed option; this directory should not be present if DWIs have already been pre-processed')
    app.console('Importing fmap data into temporary directory')
    fmap_image_list = glob.glob(os.path.join(fmap_dir, label + '*_dir-*_epi.nii*'))
    for entry in fmap_image_list:
      prefix = entry.split(os.extsep)[0]
      json_path = prefix + '.json'
      with open(json_path, 'r') as f:
        json_elements = json.load(f)
      if 'IntendedFor' in json_elements and not any(i.endswith(json_elements['IntendedFor']) for i in dwi_image_list):
        app.console('Image \'' + entry + '\' is not intended for use with DWIs; skipping')
        continue
      if not os.path.isfile(json_path):
        app.error('No sidecar JSON file found for image \'' + entry + '\'')
      # fmap files will not come with any gradient encoding in the JSON;
      #   therefore we need to add it manually ourselves so that mrcat / mrconvert can
      #   appropriately handle the table once these images are concatenated with the DWIs
      fmap_image_size = image.Header(entry).size()
      fmap_image_num_volumes = 1 if len(fmap_image_size) == 3 else fmap_image_size[3]
      fmap_dwscheme_file = 'fmap' + str(fmap_index) + '.b'
      with open(path.toTemp(fmap_dwscheme_file, False), 'w') as f:
        for line in range(0, fmap_image_num_volumes):
          f.write('0,0,1,0\n')
      run.command('mrconvert ' + entry + ' -json_import ' + json_path + ' ' + \
                  '-grad ' + path.toTemp(fmap_dwscheme_file, True) + ' ' + \
                  path.toTemp('fmap' + str(fmap_index) + '.mif', True))
      file.delTemporary(fmap_dwscheme_file)
      fmap_index += 1

    fmap_image_list = [ 'fmap' + str(index) + '.mif' for index in range(1, fmap_index) ]
  # If there's no data in fmap/ directory, need to check to see if there's any phase-encoding
  #   contrast within the input DWI(s)
  elif len(dwi_image_list) < 2 and not app.args.preprocessed:
    app.console('Inadequate data for pre-processing of subject \'' +
                label + '\': No phase-encoding contrast in input DWIs, and no fmap/ directory, '
                        'so EPI distortion correction can only be performed using Synb0-DisCo')

  dwi_image_list = [ 'dwi' + str(index) + '.mif' for index in range(1, dwi_index) ]

  # Import anatomical image
  app.console('Importing T1 image into temporary directory')
  t1w_image_list = glob.glob(os.path.join(bids_dir, label, 'anat', label + '*_T1w.nii*'))
  if len(t1w_image_list) > 1:
    app.error('More than one T1-weighted image found for subject ' + label + '; script not yet compatible with this')
  elif not t1w_image_list:
    app.error('No T1-weighted image found for subject ' + label)
  run.command('mrconvert ' + t1w_image_list[0] + ' ' + path.toTemp('T1.mif', True))

  cwd = os.getcwd()
  app.gotoTempDir()

  dwipreproc_se_epi = ''
  dwipreproc_se_epi_option = ''


  if app.args.preprocessed:

    if len(dwi_image_list) > 1:
      app.error('If DWIs have been pre-processed, then only a single DWI file should be present')
    app.console('Skipping MP-PCA denoising, Gibbs ringing removal, distortion correction and bias field correction due to use of ' + option_prefix + 'preprocessed option')
    run.function(os.rename, dwi_image_list[0], 'dwi.mif')

  else: # Do initial image pre-processing (denoising, Gibbs ringing removal if available, distortion correction & bias field correction) as normal

    # Concatenate any SE EPI images with the DWIs before denoising (& unringing), then
    #   separate them again after the fact

    if fmap_image_list:
      app.console('Concatenating DWI and fmap data for combined pre-processing')
      dwidenoise_input = 'dwi_fmap_cat.mif'
      if len(fmap_image_list) > 1:
        run.command('mrcat ' + ' '.join(fmap_image_list) + ' fmap_cat.mif -axis 3')
        fmap_header = image.Header('fmap_cat.mif')
        fmap_num_volumes = fmap_header.size()[3]
        run.command('mrcat fmap_cat.mif ' + ' '.join(dwi_image_list) + ' ' + dwidenoise_input + ' -axis 3')
        file.delTemporary('fmap_cat.mif')
      else:
        fmap_header = image.Header(fmap_image_list[0])
        fmap_num_volumes = 1 if len(fmap_header.size()) < 4 else fmap_header.size()[3]
        run.command('mrcat ' + fmap_image_list[0] + ' ' + ' '.join(dwi_image_list) + ' ' + dwidenoise_input + ' -axis 3')
      for i in fmap_image_list:
        file.delTemporary(i)
      for i in dwi_image_list:
        file.delTemporary(i)
    else:
      fmap_num_volumes = 0
      # Even if no explicit fmap images, may still need to concatenate multiple DWI inputs
      if len(dwi_image_list) > 1:
        app.console('Concatenating input DWI series')
        dwidenoise_input = 'dwi_cat.mif'
        run.command('mrcat ' + ' '.join(dwi_image_list) + ' ' + dwidenoise_input + ' -axis 3')
        for i in dwi_image_list:
          file.delTemporary(i)
      else:
        dwidenoise_input = dwi_image_list[0]


    # Step 1: Denoise
    app.console('Performing MP-PCA denoising of DWI' + (' and fmap' if fmap_num_volumes else '') + ' data')
    mrdegibbs_input = os.path.splitext(dwidenoise_input)[0] + '_denoised.mif'
    run.command('dwidenoise ' + dwidenoise_input + ' ' + mrdegibbs_input)
    file.delTemporary(dwidenoise_input)

    # Step 2: Gibbs ringing removal
    app.console('Performing Gibbs ringing removal for DWI' + (' and fmap' if fmap_num_volumes else '') + ' data')
    mrdegibbs_output = os.path.splitext(mrdegibbs_input)[0] + '_degibbs.mif'
    run.command('mrdegibbs ' + mrdegibbs_input + ' ' + mrdegibbs_output + ' -nshifts 50')
    file.delTemporary(mrdegibbs_input)

    # If fmap images and DWIs have been concatenated, now is the time to split them back apart
    if fmap_num_volumes:
      app.console('Separating DWIs and fmap images from concatenated series')
      dwipreproc_input = 'dwipreproc_in.mif'
      dwipreproc_se_epi = 'se_epi.mif'
      run.command('mrconvert ' + mrdegibbs_output + ' ' + dwipreproc_se_epi + ' -coord 3 0:' + str(fmap_num_volumes-1))
      cat_num_volumes = image.Header(mrdegibbs_output).size()[3]
      app.var(cat_num_volumes)
      run.command('mrconvert ' + mrdegibbs_output + ' ' + dwipreproc_input + ' -coord 3 ' + str(fmap_num_volumes) + ':' + str(cat_num_volumes-1))
      file.delTemporary(mrdegibbs_output)
      dwipreproc_se_epi_option = ' -se_epi ' + dwipreproc_se_epi + ' -align_seepi'
    else:
      dwipreproc_input = mrdegibbs_output
      dwipreproc_se_epi = None
      dwipreproc_se_epi_option = ''

    # Step 3: Distortion correction
    if not app.args.no_EddyandBias:
      app.console('Performing various geometric corrections of DWIs')
      dwipreproc_input_header = image.Header(dwipreproc_input)
      have_slice_timing = 'SliceTiming' in dwipreproc_input_header.keyval()
      app.var(have_slice_timing)
      mb_factor = int(dwipreproc_input_header.keyval().get('MultibandAccelerationFactor', '1'))
      app.var(mb_factor)
      if 'SliceDirection' in dwipreproc_input_header.keyval():
        slice_direction_code = dwipreproc_input_header.keyval()['SliceDirection']
        if 'i' in slice_direction_code:
          num_slices = dwipreproc_input_header.size()[0]
        elif 'j' in slice_direction_code:
          num_slices = dwipreproc_input_header.size()[1]
        elif 'k' in slice_direction_code:
          num_slices = dwipreproc_input_header.size()[2]
        else:
          num_slices = dwipreproc_input_header.size()[2]
          app.warn('Error reading BIDS field \'SliceDirection\' (value: \'' + slice_direction_code + '\'); assuming third axis')
      else:
        num_slices = dwipreproc_input_header.size()[2]
      app.var(num_slices)
      mporder = 1 + num_slices/(mb_factor*4)
      app.var(mporder)
      eddy_binary = fsl.eddyBinary(True)
      if eddy_binary:
        eddy_cuda = True
      else:
        eddy_binary = fsl.eddyBinary(False)
        eddy_cuda = False

      # modified by chenfei: comment below
      app.var(eddy_binary, eddy_cuda)
      eddy_options = []
      if have_slice_timing and eddy_cuda:
        eddy_options.append('--repol')
        eddy_options.append('--mporder=' + str(mporder))
      dwipreproc_eddy_option = ' -eddy_options \" ' + ' '.join(eddy_options) + '\"' if eddy_options else ''
      # dwipreproc_eddy_option = ''
      run.function(os.makedirs, 'eddyqc')
      # run.command('mrconvert ' + dwipreproc_input + ' dwi_preprocessed.mif') # remove dwifslpreproc to save time, modified by Chenfei
      # run.command('dwifslpreproc ' + dwipreproc_input + ' dwi_preprocessed.mif -rpe_none -pe_dir ap ' + dwipreproc_se_epi_option + dwipreproc_eddy_option + ' -eddyqc_text eddyqc/')
      if app.args.synb0: # prepare input files for synb0: create T1.nii.gz and b0.nii.gz
        run.command('mrconvert T1.mif T1.nii.gz -strides +1,+2,+3')
        run.command('dwiextract ' + dwipreproc_input + ' -bzero - | '
                      'mrcalc - 0.0 -max - | '
                      'mrmath - mean -axis 3 b0.nii.gz')
        temp_dir = os.getcwd()

        run.command('mrconvert ' + dwipreproc_input + ' dwi_denoised_degibbs.nii.gz -strides +1,+2,+3,+4 -export_grad_fsl bvecs.txt bvals.txt -json_export dwi_denoised_degibbs.json -export_pe_eddy acqparams.txt index.txt')
        # change acqparams.txt according to https://github.com/MASILab/Synb0-DISCO
        with open('acqparams.txt', 'r') as f:
          line1 = f.readline()
        with open('acqparams.txt', 'a') as f:
          f.write(line1[:6] + ' 0')
        # shutil.copytree('/home/amax/mrtrix3_connectome/synb0', os.path.join(temp_dir, 'synb0'))
        if not os.path.exists('synb0'):
          os.mkdir('synb0')
        run.command('docker run --rm -v ' + temp_dir + '/:/INPUTS/ -v ' + temp_dir + '/synb0:/OUTPUTS/ -v /usr/local/freesurfer/license.txt:/extra/freesurfer/license.txt hansencb/synb0:latest')
        file.delTemporary('T1.nii.gz')
        file.delTemporary('b0.nii.gz')

        # prepare input files for eddy: create dwi_mask
        run.command('dwi2mask ' + dwipreproc_input + ' dwi_denoised_degibbs_mask.nii.gz')
        run.command('eddy_cuda9.1 --imain=dwi_denoised_degibbs.nii.gz --mask=dwi_denoised_degibbs_mask.nii.gz --acqp=acqparams.txt --index=index.txt \
    --bvecs=bvecs.txt --bvals=bvals.txt --topup=synb0/topup --out=dwi_preprocessed ' + ' '.join(eddy_options))
        run.command('mrconvert dwi_preprocessed.nii.gz dwi_preprocessed.mif -fslgrad bvecs.txt bvals.txt -json_import dwi_denoised_degibbs.json')
      else:
        run.command('dwifslpreproc ' + dwipreproc_input + ' dwi_preprocessed.mif -rpe_header ' + dwipreproc_se_epi_option + dwipreproc_eddy_option + ' -eddyqc_text eddyqc/')

      file.delTemporary(dwipreproc_input)
      if dwipreproc_se_epi:
        file.delTemporary(dwipreproc_se_epi)

      # Step 4: Bias field correction
      if dwibiascorrect_algo:
        app.console('Performing initial B1 bias field correction of DWIs')
        # run.command('dwibiascorrect ' + dwibiascorrect_algo + ' ' + dwipreproc_input + ' dwi.mif ')  # modified by chenfei
        run.command('dwibiascorrect ' + dwibiascorrect_algo + ' dwi_preprocessed.mif dwi.mif ') # modified by chenfei
        file.delTemporary('dwi_preprocessed.mif')
      else:
        run.function(shutil.move, 'dwi_preprocessed.mif', 'dwi.mif')

    else: # no_EddyandBias
      run.command('mrconvert ' + dwipreproc_input + ' dwi.mif')



  # Step 5: Generate a brain mask for DWI
  #   Also produce a dilated version of the mask for later use
  app.console('Estimating a brain mask for DWIs')
  # # to avoid error induced by 'dwi2mask dwi.mif dwi_mask.mif' (very high signal in boundary of FOV due to synb0), create an intial mask for later use. modified by chenfei
  # run.command('maskfilter dwi_denoised_degibbs_mask.nii.gz dilate dwi_denoised_degibbs_mask_dilated.mif -npass 6')
  # run.command('mrcalc dwi.mif dwi_denoised_degibbs_mask_dilated.mif -mult dwi.mif -force')

  run.command('dwi2mask dwi.mif dwi_mask.mif')
  run.command('maskfilter dwi_mask.mif dilate dwi_mask_dilated.mif -npass 3')

  # Step 6: Generate an FA image
  # This is required for group-level analysis
  app.console('Generating FA image for group-level analysis')
  run.command('dwi2tensor dwi.mif -mask dwi_mask.mif - | tensor2metric - -fa fa.mif')

  # Step 7: Estimate response functions for spherical deconvolution
  app.console('Estimating tissue response functions for spherical deconvolution')
  run.command('dwi2response dhollander dwi.mif response_wm.txt response_gm.txt response_csf.txt -mask dwi_mask.mif')

  # Determine whether we are working with single-shell or multi-shell data
  bvalues = [int(round(float(value))) for value in image.mrinfo('dwi.mif', 'shell_bvalues').strip().split()]
  multishell = (len(bvalues) > 2)

  # Step 8: Perform spherical deconvolution
  #          Use a dilated mask for spherical deconvolution as a 'safety margin' -
  #          ACT should be responsible for stopping streamlines before they reach the edge of the DWI mask
  app.console('Estimating' + ('' if multishell else ' white matter') + ' Fibre Orientation Distribution' + ('s' if multishell else ''))
  if multishell:
    run.command('dwi2fod msmt_csd dwi.mif response_wm.txt FOD_WM.mif response_gm.txt FOD_GM.mif response_csf.txt FOD_CSF.mif '
                '-mask dwi_mask_dilated.mif -lmax 10,0,0')
    run.command('mrconvert FOD_WM.mif - -coord 3 0 | mrcat FOD_CSF.mif FOD_GM.mif - tissues.mif -axis 3')
  else:
    # Still use the msmt_csd algorithm with single-shell data: Use hard non-negativity constraint
    # Also incorporate the CSF response to provide some fluid attenuation
    run.command('dwi2fod msmt_csd dwi.mif response_wm.txt FOD_WM.mif response_csf.txt FOD_CSF.mif '
                '-mask dwi_mask_dilated.mif -lmax 10,0')
    file.delTemporary('FOD_CSF.mif')


  # Step 9: Perform brain extraction and bias field correction on the T1 image
  #         in its original space (this is necessary for histogram matching
  #         prior to registration)
  app.console('Performing brain extraction and B1 bias field correction of T1 image')
  T1_header = image.Header('T1.mif')
  T1_revert_strides_option = ' -strides ' + ','.join([str(i) for i in T1_header.strides()])
  if brain_extraction_cmd == 'runROBEX.sh':
    # Do a semi-iterative approach here: Get an initial brain mask, use that
    #   mask to estimate a bias field, then re-compute the brain mask
    run.command('mrconvert T1.mif T1.nii -strides +1,+2,+3')
    run.command(brain_extraction_cmd + ' T1.nii T1_initial_brain.nii T1_initial_mask.nii')
    file.delTemporary('T1_initial_brain.nii')
    run.command('N4BiasFieldCorrection -i T1.nii -w T1_initial_mask.nii -o T1_biascorr.nii')
    file.delTemporary('T1.nii')
    file.delTemporary('T1_initial_mask.nii')
    run.command(brain_extraction_cmd + ' T1_biascorr.nii T1_biascorr_brain.nii T1_biascorr_brain_mask.nii')
    file.delTemporary('T1_biascorr.nii')
    run.command('mrconvert T1_biascorr_brain.nii T1_biascorr_brain.mif' + T1_revert_strides_option)
    file.delTemporary('T1_biascorr_brain.nii')
    run.command('mrconvert T1_biascorr_brain_mask.nii T1_mask.mif -datatype bit' + T1_revert_strides_option)
    file.delTemporary('T1_biascorr_brain_mask.nii')
  else:
    run.command('mrconvert T1.mif T1.nii -strides -1,+2,+3')
    run.command(brain_extraction_cmd + ' -i T1.nii --noseg --nosubcortseg')
    file.delTemporary('T1.nii')
    run.command('mrconvert ' + fsl.findImage('T1.anat' + os.sep + 'T1_biascorr_brain') + ' T1_biascorr_brain.mif' + T1_revert_strides_option)
    run.command('mrconvert ' + fsl.findImage('T1.anat' + os.sep + 'T1_biascorr_brain_mask') + ' T1_mask.mif -datatype bit' + T1_revert_strides_option)
    file.delTemporary('T1.anat')

  # Step 10: Generate target images for T1->DWI registration
  app.console('Generating contrast-matched images for inter-modal registration between DWIs and T1')
  run.command('dwiextract dwi.mif -bzero - | '
              'mrcalc - 0.0 -max - | '
              'mrmath - mean -axis 3 dwi_meanbzero.mif')
  run.command('mrcalc 1 dwi_meanbzero.mif -div dwi_mask.mif -mult - | '
              'mrhistmatch nonlinear - T1_biascorr_brain.mif dwi_pseudoT1.mif -mask_input dwi_mask.mif -mask_target T1_mask.mif')
  run.command('mrcalc 1 T1_biascorr_brain.mif -div T1_mask.mif -mult - | '
              'mrhistmatch nonlinear - dwi_meanbzero.mif T1_pseudobzero.mif -mask_input T1_mask.mif -mask_target dwi_mask.mif')

  # Step 11: Perform T1->DWI registration
  #          Note that two registrations are performed: Even though we have a symmetric registration,
  #          generation of the two histogram-matched images means that you will get slightly different
  #          answers depending on which synthesized image & original image you use
  app.console('Performing registration between DWIs and T1')
  run.command('mrregister T1_biascorr_brain.mif dwi_pseudoT1.mif -type rigid -mask1 T1_mask.mif -mask2 dwi_mask.mif -rigid rigid_T1_to_pseudoT1.txt')
  file.delTemporary('T1_biascorr_brain.mif')
  run.command('mrregister T1_pseudobzero.mif dwi_meanbzero.mif -type rigid -mask1 T1_mask.mif -mask2 dwi_mask.mif -rigid rigid_pseudobzero_to_bzero.txt')
  run.command('transformcalc rigid_T1_to_pseudoT1.txt rigid_pseudobzero_to_bzero.txt average rigid_T1_to_dwi.txt')
  file.delTemporary('rigid_T1_to_pseudoT1.txt')
  file.delTemporary('rigid_pseudobzero_to_bzero.txt')
  run.command('mrtransform T1.mif T1_registered.mif -linear rigid_T1_to_dwi.txt')
  file.delTemporary('T1.mif')
  # Note: Since we're using a mask from fsl_anat (which crops the FoV), but using it as input to 5ttge fsl
  #   (which is receiving the raw T1), we need to resample in order to have the same dimensions between these two
  run.command('mrtransform T1_mask.mif T1_mask_registered.mif -linear rigid_T1_to_dwi.txt -template T1_registered.mif -interp nearest -datatype bit')
  file.delTemporary('T1_mask.mif')

  # Step 12a: Generate 5TT image for ACT
  app.console('Generating five-tissue-type (5TT) image for Anatomically-Constrained Tractography (ACT)')
  run.command('5ttgen fsl T1_registered.mif 5TT.mif -mask T1_mask_registered.mif -nocrop')
  if app.args.output_verbosity > 1:
    run.command('5tt2vis 5TT.mif vis.mif')

  # Step 12b (optional): include lesion label in 5TT image #modified by chenfei
  ## find lesion label for T2w if possible
  app.console('find lesion label for T2w if possible ')
  t2w_lesion_list = glob.glob(os.path.join(bids_dir, 'derivatives', 'lesion', label,  label + '*_T2w_label-lesion_roi.nii*'))
  if len(t2w_lesion_list) > 1:
    app.error('More than one T2-weighted lesion label found for subject ' + label + '; script not yet compatible with this')
  elif not t2w_lesion_list:
    app.console('No T2-weighted lesion label found for subject ' + label)
  else:
    # run.command('fslorient -copyqform2sform ' + t2w_lesion_list[0]) # incase qform and sform inconsistent (error would be induced in robex)
    run.command('mrconvert ' + t2w_lesion_list[0] + ' ' + path.toTemp('lesion.mif', True))
    ## Perform T2w->DWI registration
    app.console('Performing registration between T2w and DWIs')
    t2w_image_list = glob.glob(os.path.join(bids_dir, label, 'anat', label + '*_T2w.nii*'))
    t2w_image_list.sort()
    if not t2w_image_list:
      app.error('No T2-weighted image found for subject ' + label)
    else:
      # run.command('mrconvert ' + t2w_image_list[0] + ' ' + path.toTemp('T2w.mif', True))
      # run.command('mrconvert T2w.mif T2w.nii -strides -1,+2,+3')
      # run.command('fslorient -copyqform2sform ' + t2w_image_list[0])  # incase qform and sform inconsistent (error would be induced in robex)
      run.command(brain_extraction_cmd + ' ' + t2w_image_list[-1] + ' T2w_brain.nii T2w_brain_mask.nii')
      run.command('mrconvert T2w_brain.nii T2w_brain.mif -strides +1,+2,+3')
      run.command('mrconvert T2w_brain_mask.nii T2w_brain_mask.mif -strides +1,+2,+3')
      # run.command('mrregister T2w_brain.mif T1_pseudobzero.mif -mask1 T2w_brain_mask.mif -mask2 T1_mask_registered.mif -type rigid -rigid rigid_T2_to_T1.txt')
      run.command('mrregister dwi_meanbzero.mif T2w_brain.mif -mask1 dwi_mask.mif -mask2 T2w_brain_mask.mif -type rigid -rigid rigid_dwi_to_T2.txt')
      run.command('transformcompose rigid_T1_to_dwi.txt rigid_dwi_to_T2.txt rigid_T1_to_T2.txt')
      run.command('transformcalc rigid_T1_to_T2.txt invert rigid_T2_to_T1.txt')
      run.command('mrtransform T2w_brain.mif T2w_t1space.mif -linear rigid_T2_to_T1.txt -template T1_pseudobzero.mif')
      run.command('mrtransform lesion.mif lesion_t1space.mif -linear rigid_T2_to_T1.txt -interp nearest -template T1_pseudobzero.mif')
      app.console('finished T2w-qc')
      ## edit 5TT image
      run.command('5ttedit 5TT.mif -path lesion_t1space.mif 5TT.mif -force')

  # Step 13: Generate the grey matter parcellation
  #          The necessary steps here will vary significantly depending on the parcellation scheme selected

  if app.args.parcellation in [ 'desikan', 'destrieux', 'hcpmmp1' ]:

    app.console('Getting grey matter parcellation in subject space using FreeSurfer')

    run.command('mrconvert T1_registered.mif T1_registered.nii -strides +1,+2,+3') # uncomment by chenfei

    # If running in a container environment, and --debug is used (resulting in the
    #   scratch directory being a mounted drive), it's possible that attempting to
    #   construct a softlink may lead to an OSError
    # As such, run a test to determine whether or not it is possible to construct
    #   a softlink within the scratch directory; if it is not possible, revert to
    #   performing deep copies of the relevant FreeSurfer template directories
    freesurfer_template_link_function = os.symlink
    try:
      freesurfer_template_link_function(freesurfer_subjects_dir, 'test_softlink')
      file.delTemporary('test_softlink')
      app.debug('Using softlinks to FreeSurfer template directories')
    except OSError:
      app.debug('Unable to create softlinks; will perform deep copies of FreeSurfer template directories')
      freesurfer_template_link_function = shutil.copytree

    # Since we're instructing recon-all to use a different subject directory, we need to
    #   construct softlinks to a number of directories provided by FreeSurfer that
    #   recon-all will expect to find in the same directory as the overridden subject path
    for subdir in [ 'fsaverage', 'lh.EC_average', 'rh.EC_average' ]:
      run.function(freesurfer_template_link_function, os.path.join(freesurfer_subjects_dir, subdir), subdir)

    # Run FreeSurfer pipeline on this subject's T1 image
    # run.command('recon-all -sd ' + app.tempDir + ' -subjid freesurfer -i T1_registered.nii') # uncomment by chenfei
    # run.command('recon-all -sd ' + app.tempDir + ' -subjid freesurfer -all' + reconall_multithread_options) # uncomment by chenfei

    # Grab the relevant parcellation image and target lookup table for conversion
    # parc_image_path = os.path.join('freesurfer', 'mri') # uncomment by chenfei
    fs_output_dir = '/home/amax/data/cye/MScohort_BIDS_clean/derivatives/freesurfer/' + subject_label
    parc_image_path = os.path.join(fs_output_dir, 'freesurfer', 'mri')
    if app.args.parcellation == 'desikan':
      parc_image_path = os.path.join(parc_image_path, 'aparc+aseg.mgz')
    elif app.args.parcellation == 'destrieux':
      parc_image_path = os.path.join(parc_image_path, 'aparc.a2009s+aseg.mgz')
    else:
      # The HCPMMP1 parcellation is not applied as part of the recon-all command;
      #   need to explicitly map it to the subject
      # This requires SUBJECTS_DIR to be set; commands don't have a corresponding -sd option like recon-all

      # run._env['SUBJECTS_DIR'] = app.tempDir # uncomment by chenfei
      run._env['SUBJECTS_DIR'] = fs_output_dir # specific to MS-study, modified by chenfei
      # parc_image_path = os.path.join(parc_image_path, 'aparc.HCPMMP1+aseg.mgz') # uncomment by chenfei
      parc_image_path = os.path.join(fs_output_dir, 'freesurfer', 'mri', 'aparc.HCPMMP1+aseg.mgz')
      for hemi in [ 'l', 'r' ]:
        # run.command('mri_surf2surf --srcsubject fsaverage --trgsubject freesurfer --hemi ' + hemi + 'h --sval-annot ' + os.path.join(freesurfer_subjects_dir, 'fsaverage', 'label', hemi + 'h.HCPMMP1.annot') + ' --tval ' + os.path.join('freesurfer', 'label', hemi + 'h.HCPMMP1.annot'))
        run.command(
          'mri_surf2surf --srcsubject fsaverage --trgsubject freesurfer --hemi ' + hemi + 'h --sval-annot ' + os.path.join(
            freesurfer_subjects_dir, 'fsaverage', 'label', hemi + 'h.HCPMMP1.annot') + ' --tval ' + os.path.join(
            fs_output_dir, 'freesurfer', 'label', hemi + 'h.HCPMMP1.annot'))

      # run.command('mri_aparc2aseg --s freesurfer --old-ribbon --annot HCPMMP1 --o ' + parc_image_path)
      run.command('mri_aparc2aseg --s freesurfer --old-ribbon --annot HCPMMP1 --o ' + parc_image_path)

    # Perform the index conversion
    run.command('labelconvert ' + parc_image_path + ' ' + parc_lut_file + ' ' + mrtrix_lut_file + ' parc_init.mif')
    # if app.cleanup:
    #   run.function(shutil.rmtree, 'freesurfer')

    # Fix the sub-cortical grey matter parcellations using FSL FIRST
    run.command('labelsgmfix parc_init.mif T1_registered.mif ' + mrtrix_lut_file + ' parc.mif')
    file.delTemporary('parc_init.mif')

  elif app.args.parcellation in [ 'aal', 'aal2', 'craddock200', 'craddock400', 'perry512' ]:

    app.console('Registering to template and transforming grey matter parcellation back to subject space')

    # Use non-dilated brain masks for performing histogram matching & linear registration
    run.command('mrhistmatch linear T1_registered.mif -mask_input T1_mask_registered.mif ' + template_image_path + ' -mask_target ' + template_mask_path + ' - | ' \
                'mrconvert - T1_registered_histmatch.nii -strides -1,+2,+3')

    assert template_registration_software
    if template_registration_software == 'ants':

      # Use ANTs SyN for registration to template
      # From Klein et al., NeuroImage 2009:
      # Warp: ANTS 3 -m PR[htargeti.nii, hsourcei.nii, 1, 2] -o houtput transformi.nii -r Gauss[2,0] -t SyN[0.5] -i 30x99x11 -use-Histogram-Matching
      # Reslice: WarpImageMultiTransform 3 hlabeled sourcei.nii houtput labelsi.nii -R htargeti.nii htransformiWarp.nii htransformiAffine.txt --use-NN
      run.command('ANTS 3 -m PR[' + template_image_path + ', T1_registered_histmatch.nii, 1, 2] -o ANTS -r Gauss[2,0] -t SyN[0.5] -i 30x99x11 --use-Histogram-Matching')
      transformed_atlas_path = 'atlas_transformed.nii'
      run.command('WarpImageMultiTransform 3 ' + parc_image_path + ' ' + transformed_atlas_path + ' -R T1_registered_histmatch.nii -i ANTSAffine.txt ANTSInverseWarp.nii --use-NN')
      file.delTemporary(glob.glob('ANTSWarp.nii*')[0])
      file.delTemporary(glob.glob('ANTSInverseWarp.nii*')[0])
      file.delTemporary('ANTSAffine.txt')

    file.delTemporary('T1_registered_histmatch.nii')

    if parc_lut_file or mrtrix_lut_file:
      assert parc_lut_file and mrtrix_lut_file
      run.command('labelconvert ' + transformed_atlas_path + ' ' + parc_lut_file + ' ' + mrtrix_lut_file + ' parc.mif')
    else: # Not all parcellations need to go through the labelconvert step; they may already be numbered incrementally from 1
      run.command('mrconvert ' + transformed_atlas_path + ' parc.mif -strides T1_registered.mif')
    file.delTemporary(transformed_atlas_path)

  elif app.args.parcellation != 'none':
    app.error('Unknown parcellation scheme requested: ' + app.args.parcellation)

  file.delTemporary('T1_registered.nii')
  if app.args.output_verbosity > 2 and mrtrix_lut_file:
    run.command('label2colour parc.mif parcRGB.mif -lut ' + mrtrix_lut_file)


  # If no parcellation is requested, it is still possible to generate a whole-brain tractogram by
  #   explicitly providing the -streamlines option
  num_streamlines = None
  if app.args.streamlines:
    num_streamlines = app.args.streamlines
  elif app.args.parcellation != 'none':
    # If not manually specified, determine the appropriate number of streamlines based on the number of nodes in the parcellation:
    #   mean edge weight of 1,000 streamlines
    num_nodes = int(image.statistic('parc.mif', 'max'))
    num_streamlines = 500 * num_nodes * (num_nodes-1)
  if num_streamlines:

    # Step 14: Generate the tractogram
    app.console('Performing whole-brain fibre-tracking')
    tractogram_filepath = 'tractogram.tck'
    run.command('tckgen FOD_WM.mif ' + tractogram_filepath + ' -act 5TT.mif -backtrack -crop_at_gmwmi -maxlength 250 -power 0.33 ' \
                '-select ' + str(num_streamlines) + ' -seed_dynamic FOD_WM.mif')

    # Step 15: Use SIFT2 to determine streamline weights
    app.console('Running the SIFT2 algorithm to assign weights to individual streamlines')
    fd_scale_gm_option = ''
    if not multishell:
      fd_scale_gm_option = ' -fd_scale_gm'
    # If SIFT2 fails, reduce number of streamlines and try again
    while num_streamlines:
      run.command('tcksift2 ' + tractogram_filepath + ' FOD_WM.mif weights.csv -act 5TT.mif -out_mu mu.txt' + fd_scale_gm_option, False)
      if os.path.isfile('weights.csv'):
        break
      app.warn('SIFT2 failed, likely due to running out of RAM; reducing number of streamlines and trying again')
      num_streamlines = int(num_streamlines // 2)
      new_tractogram_filepath = 'tractogram_' + str(num_streamlines) + '.tck'
      run.command('tckedit ' + tractogram_filepath + ' ' + new_tractogram_filepath + ' -number ' + str(num_streamlines))
      file.delTemporary(tractogram_filepath)
      tractogram_filepath = new_tractogram_filepath
    if not num_streamlines:
      app.error('Unable to run SIFT2 algorithm for any number of streamlines')


    if app.args.output_verbosity > 2:
      # Generate TDIs:
      # - A TDI at DWI native resolution, with SIFT mu scaling, and precise mapping
      #   (for comparison to WM ODF l=0 term, to verify that SIFT2 has worked correctly)
      app.console('Producing Track Density Images (TDIs)')
      with open('mu.txt', 'r') as f:
        mu = float(f.read())
      run.command('tckmap ' + tractogram_filepath + ' -tck_weights_in weights.csv -template FOD_WM.mif -precise - | ' \
                'mrcalc - ' + str(mu) + ' -mult tdi_native.mif')
      # - Conventional TDI at super-resolution (mostly just because we can)
      run.command('tckmap ' + tractogram_filepath + ' -tck_weights_in weights.csv -template vis.mif -vox ' + ','.join([str(value/3.0) for value in image.Header('vis.mif').spacing() ]) + ' -datatype uint16 tdi_highres.mif')


  if app.args.parcellation != 'none':
    # Step 16: Generate the connectome
    #          Also get the mean length for each edge; this is the most likely alternative contrast to be useful
    app.console('Combining whole-brain tractogram with grey matter parcellation to produce the connectome')
    run.command('tck2connectome ' + tractogram_filepath + ' parc.mif connectome.csv -tck_weights_in weights.csv -symmetric -zero_diagonal -out_assignments assignments.csv')
    run.command('tck2connectome ' + tractogram_filepath + ' parc.mif meanlength.csv -tck_weights_in weights.csv -scale_length -symmetric -zero_diagonal -stat_edge mean')
    run.command('tck2connectome ' + tractogram_filepath + ' parc.mif invlength.csv -tck_weights_in weights.csv -scale_invlength  -symmetric -zero_diagonal -stat_edge sum')
    run.command('tck2connectome ' + tractogram_filepath + ' parc.mif invnodevol.csv -tck_weights_in weights.csv -scale_invnodevol -symmetric -zero_diagonal -stat_edge sum')

    if app.args.output_verbosity > 2:
      # Produce additional data that can be used for visualisation within mrview's connectome toolbar
      app.console('Generating geometric data for enhanced connectome visualisation')
      run.command('connectome2tck ' + tractogram_filepath + ' assignments.csv exemplars.tck -tck_weights_in weights.csv -exemplars parc.mif -files single')
      run.command('label2mesh parc.mif nodes.obj')
      run.command('meshfilter nodes.obj smooth nodes_smooth.obj')
      file.delTemporary('nodes.obj')


  # Prepare output path for writing
  app.console('Processing for subject \'' + label + '\' completed; writing results to output directory')
  if os.path.exists(output_dir):
    run.function(shutil.rmtree, output_dir)
  run.function(os.makedirs, output_dir)
  run.function(os.makedirs, os.path.join(output_dir, 'connectome'))
  run.function(os.makedirs, os.path.join(output_dir, 'dwi'))
  run.function(os.makedirs, os.path.join(output_dir, 'tractogram'))
  if app.args.output_verbosity > 1:
    run.function(os.makedirs, os.path.join(output_dir, 'anat'))

  parc_string = '_parc-' + app.args.parcellation

  # Generate a copy of the lookup table file:
  #   - Use the post-labelconvert file if it's used; otherwise, if the atlas
  #     itself comes with a lookup table that didn't require conversion, write that;
  #   - In the group directory rather than the subject directory;
  #   - If it doesn't already exist.
  lut_export_file = mrtrix_lut_file if mrtrix_lut_file else parc_lut_file
  if lut_export_file:
    lut_export_path = os.path.join(output_prefix, parc_string[1:] + '_lookup' + os.path.splitext(lut_export_file)[1])
    try:
      shutil.copy(lut_export_file, lut_export_path)
    except OSError:
      pass

  # Copy / convert necessary files to output directory
  if app.args.parcellation != 'none':
    run.function(shutil.copy, 'connectome.csv', os.path.join(output_dir, 'connectome', label + parc_string + '_level-participant_connectome.csv'))
  if num_streamlines:
    run.function(shutil.copy, 'mu.txt', os.path.join(output_dir, 'tractogram', label + '_mu.txt'))
  run.command('mrconvert dwi_meanbzero.mif ' + os.path.join(output_dir, 'dwi', label + '_meanbzero.nii.gz') + ' -strides +1,+2,+3')
  run.command('mrconvert fa.mif ' + os.path.join(output_dir, 'dwi', label + '_model-tensor_fa.nii.gz') + ' -strides +1,+2,+3')
  run.function(shutil.copy, 'response_wm.txt', os.path.join(output_dir, 'dwi', label + '_tissue-WM_response.txt'))
  with open(os.path.join(output_dir, 'dwi', label + '_bvalues.txt'), 'w') as f:
    f.write(' '.join([str(value) for value in bvalues]))
  if app.args.output_verbosity > 1:
    if app.args.parcellation != 'none':
      run.command('mrconvert parc.mif ' + os.path.join(output_dir, 'anat', label + parc_string + '_indices.nii.gz') + ' -strides +1,+2,+3')
      run.function(shutil.copy, 'meanlength.csv', os.path.join(output_dir, 'connectome', label + parc_string + '_meanlength.csv'))
      run.function(shutil.copy, 'invlength.csv', os.path.join(output_dir, 'connectome', label + parc_string + '_invlength.csv'))
      run.function(shutil.copy, 'invnodevol.csv', os.path.join(output_dir, 'connectome', label + parc_string + '_invnodevol.csv'))
    run.command('mrconvert dwi.mif ' + os.path.join(output_dir, 'dwi', label + '_dwi.nii.gz') + \
                ' -export_grad_fsl ' + os.path.join(output_dir, 'dwi', label + '_dwi.bvec') + ' ' + os.path.join(output_dir, 'dwi', label + '_dwi.bval') + \
                ' -strides +1,+2,+3,+4')
    run.command('mrconvert dwi_mask.mif ' + os.path.join(output_dir, 'dwi', label + '_brainmask.nii.gz') + ' -datatype uint8 -strides +1,+2,+3')
    run.command('mrconvert T1_registered.mif ' + os.path.join(output_dir, 'anat', label + '_T1w.nii.gz') + ' -strides +1,+2,+3')
    run.command('mrconvert T1_mask_registered.mif ' + os.path.join(output_dir, 'anat', label + '_brainmask.nii.gz') + ' -datatype uint8 -strides +1,+2,+3')
    run.command('mrconvert 5TT.mif ' + os.path.join(output_dir, 'anat', label + '_5TT.nii.gz') + ' -strides +1,+2,+3,+4')
    run.command('mrconvert vis.mif ' + os.path.join(output_dir, 'anat', label + '_tissues3D.nii.gz') + ' -strides +1,+2,+3')
    run.command('mrconvert FOD_WM.mif ' + os.path.join(output_dir, 'dwi', label + '_tissue-WM_ODF.nii.gz') + ' -strides +1,+2,+3,+4')
    if multishell:
      run.function(shutil.copy, 'response_gm.txt', os.path.join(output_dir, 'dwi', label + '_tissue-GM_response.txt'))
      run.function(shutil.copy, 'response_csf.txt', os.path.join(output_dir, 'dwi', label + '_tissue-CSF_response.txt'))
      run.command('mrconvert FOD_GM.mif ' + os.path.join(output_dir, 'dwi', label + '_tissue-GM_ODF.nii.gz') + ' -strides +1,+2,+3,+4')
      run.command('mrconvert FOD_CSF.mif ' + os.path.join(output_dir, 'dwi', label + '_tissue-CSF_ODF.nii.gz') + ' -strides +1,+2,+3,+4')
      run.command('mrconvert tissues.mif ' + os.path.join(output_dir, 'dwi', label + '_tissue-all.nii.gz') + ' -strides +1,+2,+3,+4')
    if not app.args.preprocessed:
      run.function(shutil.copytree, 'eddyqc', os.path.join(output_dir, 'dwi', 'eddyqc'))
  if app.args.output_verbosity > 2:
    if t2w_lesion_list:
      run.command('mrconvert T2w_t1space.mif ' + os.path.join(output_dir, 'anat', label + '_T2w_t1space.nii.gz') + ' -strides +1,+2,+3')
      run.command('mrconvert T2w_brain.mif ' + os.path.join(output_dir, 'anat', label + '_T2w_brain.nii.gz') + ' -strides +1,+2,+3')
      run.command('mrconvert lesion.mif ' + os.path.join(output_dir, 'anat', label + '_lesion.nii.gz') + ' -strides +1,+2,+3')
      run.command('mrconvert lesion_t1space.mif ' + os.path.join(output_dir, 'anat', label + '_lesion_t1space.nii.gz') + ' -strides +1,+2,+3')
    if num_streamlines:
      # Move rather than copying the tractogram just because of its size
      run.function(shutil.move, tractogram_filepath, os.path.join(output_dir, 'tractogram', label + '_tractogram.tck'))
      run.function(shutil.copy, 'weights.csv', os.path.join(output_dir, 'tractogram', label + '_weights.csv'))
      run.command('mrconvert tdi_native.mif ' + os.path.join(output_dir, 'tractogram', label + '_variant-native_tdi.nii.gz') + ' -strides +1,+2,+3')
      run.command('mrconvert tdi_highres.mif ' + os.path.join(output_dir, 'tractogram', label + '_variant-highres_tdi.nii.gz') + ' -strides +1,+2,+3')
    if app.args.parcellation != 'none':
      run.function(shutil.copy, 'assignments.csv', os.path.join(output_dir, 'connectome', label + parc_string + '_assignments.csv'))
      run.function(shutil.copy, 'exemplars.tck', os.path.join(output_dir, 'connectome', label + parc_string + '_exemplars.tck'))
      run.function(shutil.copy, 'nodes_smooth.obj', os.path.join(output_dir, 'anat', label + parc_string + '.obj'))
      run.command('mrconvert parcRGB.mif ' + os.path.join(output_dir, 'anat', label + parc_string +'_colour.nii.gz') + ' -strides +1,+2,+3')

  # Manually wipe and zero the temp directory (since we might be processing more than one subject)
  os.chdir(cwd)
  if app.cleanup:
    app.console('Deleting temporary directory ' + app.tempDir)
    # Can't use run.function() here; it'll try to write to the log file that resides in the temp directory just deleted
    shutil.rmtree(app.tempDir)
  else:
    app.console('Contents of temporary directory kept; location: ' + app.tempDir)
  app.tempDir = ''

# End of runSubject() function


analysis_choices = [ 'participant', 'group' ]
parcellation_choices = [ 'aal', 'aal2', 'craddock200', 'craddock400', 'desikan', 'destrieux', 'hcpmmp1', 'none', 'perry512' ]
registration_choices = [ 'ants', 'fsl' ]

app.init('cye', 'Generate structural connectomes based on diffusion-weighted and T1-weighted image data using state-of-the-art reconstruction tools, particularly those provided in MRtrix3')
app.cmdline._action_groups[2].add_argument(option_prefix + 'skip-bids-validator', dest='skipbidsvalidator', action='store_false', help='Skip BIDS validation')
app.cmdline.add_argument('bids_dir', help='The directory with the input dataset formatted according to the BIDS standard.')
app.cmdline.add_argument('output_dir', help='The directory where the output files should be stored. If you are running group level analysis, this folder should be prepopulated with the results of the participant level analysis.')
app.cmdline.add_argument('analysis_level', help='Level of the analysis that will be performed. Multiple participant level analyses can be run independently (in parallel) using the same output_dir. Options are: ' + ', '.join(analysis_choices), choices=analysis_choices)
batch_options = app.cmdline.add_argument_group('Options specific to the batch processing of subject data')
batch_options.add_argument(option_prefix + 'participant_label', nargs='+', help='The label(s) of the participant(s) that should be analyzed.'
                                                                                ' The label(s) correspond(s) to sub-<participant_label> from the BIDS spec (so it does _not_ include "sub-"). '
                                                                                'If this parameter is not provided, all subjects will be analyzed sequentially. Multiple participants can be specified with a space-separated list.')
participant_options = app.cmdline.add_argument_group('Options that are relevant to participant-level analysis')

participant_options.add_argument(option_prefix + 'atlas_path', metavar='path', help='The filesystem path in which to search for an atlas parcellation (may be necessary if not installed in the same location as they are placed within the BIDS App container)')
participant_options.add_argument(option_prefix + 'output_verbosity', type=int, default=3, help='The verbosity of script output (number from 1 to 3); higher values result in more generated data being included in the output directory')

participant_options.add_argument(option_prefix + 'parcellation', help='The choice of connectome parcellation scheme (compulsory for participant-level analysis). Options are: ' + ', '.join(parcellation_choices), choices=parcellation_choices)
participant_options.add_argument(option_prefix + 'preprocessed', action='store_true', help='Indicate that the subject DWI data have been preprocessed, and hence initial image processing steps will be skipped')
participant_options.add_argument(option_prefix + 'streamlines', type=int, help='The number of streamlines to generate for each subject (will be determined heuristically if not explicitly set)')
participant_options.add_argument(option_prefix + 'template_reg', metavar='software', help='The choice of registration software for mapping subject to template space. Options are: ' + ', '.join(registration_choices), choices=registration_choices)
participant_options.add_argument(option_prefix + 'synb0', action='store_true', help='Indicate that synb0 docker image was used to simulate fsl-topup')
participant_options.add_argument(option_prefix + 'no_EddyandBias', action='store_true', help='Indicate that eddy-correction/topup/dwibiascorrect were not used, just for quick-test')


app.parse()


if app.isWindows():
  app.error('Script cannot be run on Windows due to FSL dependency')

if app.args.skipbidsvalidator:
  app.console('Skipping BIDS validation based on user request')
elif find_executable('bids-validator'):
  run.command('bids-validator ' + app.args.bids_dir)
else:
  app.warn('BIDS validator script not installed; proceeding without validation of input data')

# Running participant level
if app.args.analysis_level == 'participant':

  if app.args.output_verbosity < 1 or app.args.output_verbosity > 3:
    app.error('Valid values for ' + option_prefix + 'output_verbosity option are from 1 to 3')

  subjects_to_analyze = [ ]
  # Only run a subset of subjects
  if app.args.participant_label:
    subjects_to_analyze = [ 'sub-' + sub_index for sub_index in app.args.participant_label ]
    for subject_dir in subjects_to_analyze:
      if not os.path.isdir(os.path.join(app.args.bids_dir, subject_dir)):
        app.error('Unable to find directory for subject: ' + subject_dir)
  # Run all subjects sequentially
  else:
    subject_dirs = glob.glob(os.path.join(app.args.bids_dir, 'sub-*'))
    subjects_to_analyze = [ 'sub-' + directory.split("-")[-1] for directory in subject_dirs ]
    if not subjects_to_analyze:
      app.error('Could not find any subjects in BIDS directory')

  for subject_label in subjects_to_analyze:
    app.console('Commencing execution for subject ' + subject_label)
    runSubject(app.args.bids_dir, subject_label, os.path.abspath(app.args.output_dir))


app.complete()
