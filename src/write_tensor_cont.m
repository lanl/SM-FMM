function write_tensor_cont (T1,T1dim,T2,T2dim,icomplex,fidout)
  
  if (length(T1) ~= (length(T1dim) + 1))
    error('length(T1) ~= length(T1dim) + 1')
  endif
  
  if (length(T2) ~= (length(T2dim) + 1))
    error('length(T2) ~= length(T2dim) + 1')
  endif
 
  indices_all = [T1(2:end),T2(2:end)];
  dim_all = [T1dim,T2dim];
  indices = unique(indices_all);
  repeated_ind = '';
  repeated_ind_dim = '';
  
  # remove repeated indices
  Toutdim = dim_all;
  Tout = indices_all;
  ind_rep = [];
  for i = 1:length(indices)
    if (sum(indices_all == indices(i)) == 2)
      
      ind2 = find(indices_all == indices(i));
      if (dim_all(ind2(1)) ~= dim_all(ind2(2)))
        error('problem in repeated indices dimension')
      endif
      ind_rep = [ind_rep,ind2];
 
      repeated_ind = [repeated_ind, indices(i)];
      repeated_ind_dim = [repeated_ind_dim,dim_all(ind2(1))];
  
    else
  
      if (sum(indices_all ~= indices(i)) == 1)
        error('problem in indices')
      endif 
  
    endif
  endfor
  Toutdim(ind_rep) = [];
  Tout(ind_rep) = [];
  Tout = ['T',Tout];
  
  # Tout
  # Toutdim
  # repeated_ind 
  # indices_all
  
  if (icomplex == 0)
    str = ['  pure function ',T1,T2,'(',T1,',',T2,') result (',Tout,')'];
    fprintf(fidout, [str, '\n']);
  elseif (icomplex == 1)
    str = ['  pure function ',T1,T2,'_cmplx(',T1,',',T2,') result (',Tout,')'];
    fprintf(fidout, [str, '\n']);
  endif
  
  str = ['    implicit none'];
  fprintf(fidout, [str, '\n']);
  
  tmp1 = T1dim(1);
  for i = 2:length(T1dim)
    tmp1 = [tmp1,',',T1dim(i)];
  endfor
  tmp2 = T2dim(1);
  for i = 2:length(T2dim)
    tmp2 = [tmp2,',',T2dim(i)];
  endfor
  if (icomplex == 0)
    str = ['    double precision, intent(in) :: ',T1,'(',tmp1,'), ',T2,'(',tmp2,')'];
  elseif (icomplex == 1)
    str = ['    double complex, intent(in) :: ',T1,'(',tmp1,'), ',T2,'(',tmp2,')'];
  endif
  fprintf(fidout, [str, '\n']);
  
  tmpout = Toutdim(1);
  for i = 2:length(Toutdim)
    tmpout = [tmpout,',',Toutdim(i)];
  endfor
  if (icomplex == 0)
    str = ['    double precision :: ',Tout,'(',tmpout,'), dum'];
  elseif (icomplex == 1)
    str = ['    double complex :: ',Tout,'(',tmpout,'), dum'];
  endif
  fprintf(fidout, [str, '\n']);
  
  tmp_ind = indices(1);
  for i = 2:length(indices)
    tmp_ind = [tmp_ind,',',indices(i)];
  endfor
  str = ['    integer :: ',tmp_ind];
  fprintf(fidout, [str, '\n']);
  
  str = '';
  fprintf(fidout, [str, '\n']);
  
  for i = 1:length(Toutdim)
  
    str = ['    do ',Tout(i+1),' = 1,',Toutdim(i)];
    fprintf(fidout, [str, '\n']);
  
  end
  
  str = '      dum = 0.0';
  fprintf(fidout, [str, '\n']);
  
  for i = 1:length(repeated_ind)
  
    str = ['      do ',repeated_ind(i),' = 1,',repeated_ind_dim(i)];
    fprintf(fidout, [str, '\n']);
  
  end
  
  tmp1 = T1(2);
  for i = 3:length(T1)
    tmp1 = [tmp1,',',T1(i)];
  endfor
  tmp2 = T2(2);
  for i = 3:length(T2)
    tmp2 = [tmp2,',',T2(i)];
  endfor
  str = ['        dum = dum + ',T1,'(',tmp1,')*',T2,'(',tmp2,')'];
  fprintf(fidout, [str, '\n']);
  
  for i = 1:length(repeated_ind)
  
    str = '      enddo';
    fprintf(fidout, [str, '\n']);
  
  end
  
  tmpout = Tout(2);
  for i = 3:length(Tout)
    tmpout = [tmpout,',',Tout(i)];
  endfor
  str = ['      ',Tout,'(',tmpout,') = dum'];
  fprintf(fidout, [str, '\n']);
  
  for i = 1:length(Toutdim)
  
    str = '    enddo';
    fprintf(fidout, [str, '\n']);
  
  end
  
  str = '';
  fprintf(fidout, [str, '\n']);
  
  str = '    return';
  fprintf(fidout, [str, '\n']);
  
  if (icomplex == 0)
    str = ['  end function ',T1,T2];
    fprintf(fidout, [str, '\n']);
  elseif (icomplex == 1)
    str = ['  end function ',T1,T2,'_cmplx'];
    fprintf(fidout, [str, '\n']);
  endif

  str = '';
  fprintf(fidout, [str, '\n']);

endfunction