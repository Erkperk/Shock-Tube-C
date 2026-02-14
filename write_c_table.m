function write_c_table(fid, name, data, nrows, ncols)
  fprintf(fid, 'static const double %s[%d][%d] = {\n', name, nrows, ncols);
  for i = 1:nrows
    fprintf(fid, '  {');
    for j = 1:ncols
      if j > 1
        fprintf(fid, ',');
      end
      fprintf(fid, '%.10g', data(i,j));
    end
    if i < nrows
      fprintf(fid, '},\n');
    else
      fprintf(fid, '}\n');
    end
  end
  fprintf(fid, '};\n\n');
end
