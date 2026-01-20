clear;
if ~exist('results', 'dir')
  mkdir('results')
end
datetime.setDefaultFormats('default','_yyyy-MM-dd_HH-mm');
current_dt = string(datetime('now'));
[filename, location] = uigetfile({'*.xlsx'},'Input Data File Selector','MultiSelect','on');
if ~isequal(filename, 0)
  wd = append('results/Date', current_dt, '/');
  if iscell(filename)
    result = cell(length(filename), 1);
    for i = 1:1:length(filename)
      io.wd = append(wd, '/', filename{i}, '/');
      io.filename = filename{i};
      disp(append('File #', string(i), ': ', filename{i}));
      result{i} = calculateFileData(filename{i}, location, io);
    end
  else
    result = cell(1, 1);
    io.wd = wd;
    io.filename = filename;
    result{1} = calculateFileData(filename, location, io);
  end
  save(append(wd, 'results.mat'), "result");
  visualizeTGResult(result);
end

function output_result = calculateFileData(filename, location, io)
  datetime.setDefaultFormats('default','yyyy-MM-dd HH:mm:ss');
  output_result.is_ok = false;
  input_data = readtable(fullfile(location, filename), 'ReadVariableNames', false, 'Sheet', 'Data');
  if height(input_data) > 0
    if width(input_data) > 21
      if isa(input_data.Var2, 'double') && isa(input_data.Var3, 'double') && isa(input_data.Var4, 'double') && ...
         isa(input_data.Var5, 'double') && isa(input_data.Var6, 'double') && isa(input_data.Var7, 'double') && ...
         isa(input_data.Var8, 'double') && isa(input_data.Var9, 'double') && isa(input_data.Var10, 'double') && ...
         isa(input_data.Var11, 'double') && isa(input_data.Var12, 'double') && isa(input_data.Var13, 'double') && ...
         isa(input_data.Var14, 'double') && isa(input_data.Var15, 'double') && isa(input_data.Var16, 'double') && ...
         isa(input_data.Var17, 'double') && isa(input_data.Var18, 'double') && isa(input_data.Var19, 'double') && ...
         isa(input_data.Var20, 'double') && isa(input_data.Var21, 'double')
        output_result.is_ok = true;
        grating = table2struct(table(string(input_data.Var1), input_data.Var2 * 1E-3, input_data.Var3, ...
                                     'VariableNames', {'slit' 'p' 'df'}));
        pump = table2struct(table(input_data.Var4 * pi, input_data.Var5 * 1E-3, input_data.Var6, ...
                                  input_data.Var7, input_data.Var8 * 1E-6, input_data.Var9, ...
                                  'VariableNames', {'pd' 'angle' 'wfc' 'waist' 'lambda' 'intensity'}));
        probe = table2struct(table(input_data.Var10 * pi, input_data.Var11 * 1E-3, input_data.Var12, ...
                                   input_data.Var13, input_data.Var14 * 1E-6, input_data.Var15, ...
                                   'VariableNames', {'pd' 'angle' 'wfc' 'waist' 'lambda' 'intensity'}));
        sample = table2struct(table(input_data.Var16, complex(input_data.Var17, input_data.Var18), ...
                                    'VariableNames', {'t' 'n2'}));
        if width(input_data) > 22 && isa(input_data.Var23, 'double')
          sys = table2struct(table(input_data.Var19, input_data.Var20 * 1E3, input_data.Var21, ...
                                   input_data.Var22 * 1E-3, input_data.Var23, ...
                                   'VariableNames', {'g2s' 'ff' 'psd' 'ds' 'dsdp'}));
        else
          sys = table2struct(table(input_data.Var19, input_data.Var20 * 1E3, input_data.Var21, ...
                                   input_data.Var22 * 1E-3, input_data.Var22 / 0, ...
                                   'VariableNames', {'g2s' 'ff' 'psd' 'ds' 'dsdp'}));
        end
        if ~exist(io.wd, 'dir')
          mkdir(io.wd)
        end
        status = copyfile(fullfile(location, filename), io.wd);
        io.filename = append(io.wd, filename);
        output_result.wd = io.wd;
        output_result.filename = io.filename;
        if status == true
          calculated_result = cell(height(input_data), 1);
          calc_string = ["Line #", "Status", "Results Tab", "Start Time", "Finish Time", "Notes"];
          writematrix(calc_string, io.filename, 'Sheet', 'Result', 'Range', 'A1');
          for i = 1:1:height(input_data)
            io.i = string(i);
            result.i = i;
            save_cell = string(2 + i);
            disp(append('  Line #', string(i), ': calculation started at ', string(datetime('now'))));
            start_time = string(datetime('now'));
            io.ddp = false;
            io.sdp = false;
            if (width(input_data) > 23 && contains(input_data.Var24(i), 'Y'))
              io.ddp = true;
            end
            if (width(input_data) > 24 && contains(input_data.Var25(i), 'Y'))
              io.sdp = true;
            end
            % parallel calculation can be started here
            [is_ok, x, intensity, message, io_out, pump_out, probe_out, sample_out] = ...
              processTGCalculation(grating(i), pump(i), probe(i), sample(i), sys(i), io);
            if ~is_ok
              result.is_ok = 'Not ';
            else
              result.is_ok = '';
            end
            result.is_ok = append(result.is_ok, 'OK');
            result.grating = grating(i);
            result.pump = pump_out;
            result.probe = probe_out;
            result.sample = sample_out;
            result.sys = sys(i);
            result.io = io_out;
            if is_ok
              result.result.x = x;
              result.result.intensity = intensity;
            else
              result.message = message;
              output_result.is_ok = false;
            end
            calculated_result{i} = result;
            finish_time = string(datetime('now'));
            if io_out.sdp
              calc_string = [io.i, result.is_ok, io_out.line, start_time, finish_time, message];
            else
              calc_string = [io.i, result.is_ok, "", start_time, finish_time, message];
            end
            try
              writematrix(calc_string, io.filename, 'Sheet', 'Result', 'Range', append('A', save_cell));
            catch
              warning('Datafile IO porblems. Line is skipped.');
            end
            disp(append('  Line #', string(i), ': calculation finished at ', string(datetime('now'))));
          end
        else
          calculated_result = 'Copy warning: filecopy is not possible. Please check access rights.';
          warning(calculated_result);
        end
      else
        calculated_result = 'Wrong data: data format in table has wrong type. Please check input data.';
        warning(calculated_result);
      end
    else
      calculated_result = 'Wrong fileformat: datafile is corrupted. Please check input data.';
      warning(calculated_result);
    end
  else
    calculated_result = 'Wrong data: datafile is emplty. Please add input data.';
    warning(calculated_result);
  end
  output_result.data = calculated_result;
end