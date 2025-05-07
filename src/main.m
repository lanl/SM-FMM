clear all; close all;

funct = {
'TijkTjk',
'TijkTj',
'TijklTjkl',
'TijklTjk',
'TijklTj',
'TijklmTjklm',
'TijklmTjlm',
'TijklmTjm',
'TijklmTj',
'TijklmnTjklmn',
'TijklmnTjlmn',
'TijklmnTjmn',
'TijklmnTjn',
'TijklmnTj',

'TijkTjk_cmplx',
'TijkTj_cmplx',
'TijklTjkl_cmplx',
'TijklTjk_cmplx',
'TijklTj_cmplx',
'TijklmTjklm_cmplx',
'TijklmTjlm_cmplx',
'TijklmTjm_cmplx',
'TijklmTj_cmplx',
'TijklmnTjklmn_cmplx',
'TijklmnTjlmn_cmplx',
'TijklmnTjmn_cmplx',
'TijklmnTjn_cmplx',
'TijklmnTj_cmplx',
}
fidout = fopen('tensor_contr_1d.f90', 'w');

for iarr = 1:length(funct)

  str = funct{iarr}
  
  if (str(end-4:end) == 'cmplx')
    icomplex = 1
    str(end-5:end) = []
  else
    icomplex = 0
  end
  
  ind = find(str == 'T');
  T1 = str(1:ind(2)-1)
  T2 = str(ind(2):end)
  T1dim = repmat('3',1,length(T1)-1)
  T2dim = repmat('3',1,length(T2)-1)

  write_tensor_cont_1d (T1,T1dim,T2,T2dim,icomplex,fidout)

end

fclose(fidout)
return

funct = {
'TijklmTklm',
'TmlkjiTmlk_cmplx',
'TijklmTkl',
'TijkTjk',
'TijklTjkl',
'TijkTj',
'TkjiTj_cmplx',
'TijklTjk',
'TlkjiTkj_cmplx',
'TijklTj',
'TlkjiTj_cmplx',
'TmlkjiTlk_cmplx',
'TijkTk',
'TijklmnTkl',
'TnmlkjiTlk_cmplx',
'TijklmnTklm',
'TnmlkjiTmlk_cmplx',
'TijklmnTkln',
'TnmlkjiTnlk_cmplx',
'TijklmnTklmn',
'TnmlkjiTnmlk_cmplx',
'TijklTkl',
'TlkjiTlk_cmplx',
'TijklTklmn',
'TijklmTjklm',
'TijklmnTjklmn',
'TijklmTjlm',
'TmlkjiTmlj_cmplx',
'TijklmTjm',
'TmlkjiTmj_cmplx',
'TijklmTj',
'TmlkjiTj_cmplx',
'TijklmnTjlmn',
'TijklmnTjmn',
'TijklmnTjn',
'TijklmnTj',
'TkjiTkj_cmplx',
'TlkjiTlkj_cmplx',
'TmlkjiTmlkj_cmplx',
'TnmlkjiTnmlkj_cmplx',
'TnmlkjiTnmlj_cmplx',
'TnmlkjiTnmj_cmplx',
'TnmlkjiTnj_cmplx'
'TnmlkjiTj_cmplx'
}

fidout = fopen('tensor_contr.f90', 'w');

for iarr = 1:length(funct)

  str = funct{iarr}
  
  if (str(end-4:end) == 'cmplx')
    icomplex = 1
    str(end-5:end) = []
  else
    icomplex = 0
  end
  
  ind = find(str == 'T');
  T1 = str(1:ind(2)-1)
  T2 = str(ind(2):end)
  T1dim = repmat('3',1,length(T1)-1)
  T2dim = repmat('3',1,length(T2)-1)

  write_tensor_cont (T1,T1dim,T2,T2dim,icomplex,fidout)

end

fclose(fidout)
return



# T1 = 'Tijko'; # test
# T1dim = '4734';
# T2 = 'Tjklo';
# T2dim = '7354';
