 eeglab
 EEG = pop_loadset('filename',eeg_file,'filepath',eeg_path);
 [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
 % set output folder (of), subject_name, session_name, elec_file, number of
 % layers (nl), plotting, and comp_index
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
 % Realistic modeling from MRI
 
 % Do segmentation using the GUI
 nft_mesh_generation(subject_name, of, nl)
 nft_source_space_generation(subject_name, of)
 % Do co-registration using the GUI
 nft_forward_problem_solution(subject_name, session_name, of);
 dip1 = nft_inverse_problem_solution(subject_name, session_name, of, EEG, comp_index, plotting, elec_file)
  
 % Set NFT dipole structure to EEGLAB dipole structure
lof = length(of); if of(lof) ~= filesep of(lof+1) = filesep; end

 mri_file = [of subject_name '_mri.mat'];
 vol_file = [of subject_name '_vol.mat'];
 EEG.dipfit.hdmfile = vol_file;
 EEG.dipfit.mrifile =  mri_file;
 EEG.dipfit.coordformat = 'MRI';
 EEG.dipfit.model = EEG.etc.nft.model;
    
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
 
 % Warping mesh
 nft_warping_mesh(subject_name, session_name, elec_file, nl, of,0,0);
 nft_forward_problem_solution(subject_name, session_name, of);
 dip1 = nft_inverse_problem_solution(subject_name, session_name, of, EEG, comp_index, plotting, elec_file)
    
 % Set NFT dipole structure to EEGLAB dipole structure
 eeglab_folder = dirname(which('eeglab'));
 mri_file = [eeglab_folder '/plugins/dipfit2.2/standard_BEM/standard_mri.mat'];
 vol_file = [eeglab_folder '/plugins/dipfit2.2/standard_BEM/standard_vol.mat'];
 EEG.dipfit.hdmfile = vol_file;
 EEG.dipfit.mrifile =  mri_file;
 EEG.dipfit.coordformat = 'MNI';
 EEG.dipfit.model = EEG.etc.nft.model;
 
