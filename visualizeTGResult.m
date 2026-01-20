function visualizeTGResult(tg_calc_result)
  % Postcalculation part of the code has to be implemented here
  visualizeDSDP(tg_calc_result);
  visualizeTG(tg_calc_result);
  visualizeFF(tg_calc_result);
  % End of postcalculation part of the code
end

function visualizeDSDP(tg_calc_result)
  for i = 1:1:length(tg_calc_result)
    tg_data_line = tg_calc_result{i}.data;
    if tg_calc_result{i}.is_ok
      if contains(tg_calc_result{i}.filename, 'sdp') && length(tg_data_line) > 4
        tg_i_err = zeros(length(tg_data_line) - 1, 1);
        tg_n_err = zeros(length(tg_data_line) - 1, 1);
        for j = 2:1:length(tg_data_line)
          tg_n_err(j - 1) = tg_data_line{j - 1}.sys.dsdp;
          tg_i_diff = abs(tg_data_line{j - 1}.result.intensity.x - tg_data_line{j}.result.intensity.x);
          tg_i_err(j - 1) = sum(tg_i_diff);
        end
        display1DError(tg_i_err, tg_n_err, 1e-5, 6);
        writematrix('Division parameter', tg_calc_result{i}.filename, 'Sheet', 'CalculationError', 'Range', 'A1');
        writematrix('Error', tg_calc_result{i}.filename, 'Sheet', 'CalculationError', 'Range', 'B1');
        writematrix(tg_n_err, tg_calc_result{i}.filename, 'Sheet', 'CalculationError', 'Range', 'A2');
        writematrix(tg_i_err, tg_calc_result{i}.filename, 'Sheet', 'CalculationError', 'Range', 'B2');
        saveas(gcf, append(tg_calc_result{i}.wd, 'dsdp_scale.png'));
      end
    end
  end
end

function visualizeTG(tg_calc_result)
  for i = 1:1:length(tg_calc_result)
    tg_data_line = tg_calc_result{i}.data;
    if tg_calc_result{i}.is_ok
      if contains(tg_calc_result{i}.filename, 'tg') && length(tg_data_line) > 2 && isfield(tg_data_line{1}.sample, 'tg')
        if length(tg_data_line{1}.sample.tg.x) < length(tg_data_line{2}.sample.tg.x)
          tg_x = tg_data_line{1}.sample.tg.x;
        else
          tg_x = tg_data_line{2}.sample.tg.x;
        end
        tg_y = zeros(length(tg_data_line), 1);
        tg_i = zeros(length(tg_x), length(tg_data_line));

        if contains(tg_calc_result{i}.filename, 'fd')
          for j = 1:1:length(tg_data_line)
            tg_y(j) = 1 / tg_data_line{j}.probe.wfc;
            if length(tg_x) > tg_data_line{j}.sample.tg.pd
              tg_i(1:length(tg_data_line{j}.sample.tg.pd), j) = tg_data_line{j}.sample.tg.pd;
            else
              tg_i(:, j) = tg_data_line{j}.sample.tg.pd(1:length(tg_x));
            end
          end
          y_label = 'Curvature radius [mm]';
        elseif contains(tg_calc_result{i}.filename, 'pp')
          for j = 1:1:length(tg_data_line)
            tg_y(j) = tg_data_line{j}.sample.tg.pp;
            tg_i(:, j) = tg_data_line{j}.sample.tg.pd(1:length(tg_x));
          end
          y_label = 'Pump to probe ratio';
        elseif contains(tg_calc_result{i}.filename, 'ib')
          for j = 1:1:length(tg_data_line)
            tg_y(j) = tg_data_line{j}.sample.tg.beta;
            tg_i(:, j) = tg_data_line{j}.result.intensity.x(1:length(tg_x));
          end
          y_label = 'Phase effectivity [rad]';
        elseif contains(tg_calc_result{i}.filename, 'g2s')
          %fid = fopen(append(tg_calc_result{i}.wd, 'result_tg_g2s.txt'), 'wt');
          for j = 1:1:length(tg_data_line)
            tg_y(j) = tg_data_line{j}.sys.g2s;
            tg_i(:, j) = tg_data_line{j}.sample.tg.pd(1:length(tg_x));
            %for k = 1:1:length(tg_x)
            %  fprintf(fid, '%f %f %f\n', tg_y(j), tg_x(k), tg_i(k, j));
            %end
          end
          %fclose(fid);
          y_label = 'Distance [mm]';
        end

        display2DImage(tg_i, tg_x, tg_y, 'TG image', 'x [{\mu}m]', y_label, 7);
        imwrite(tg_i * 255, jet, append(tg_calc_result{i}.wd, 'tg_2d.png'), 'png');
        saveas(gcf, append(tg_calc_result{i}.wd, 'tg_2d_scale.png'));
      end
    end
  end
end

function visualizeFF(tg_calc_result)
  for i = 1:1:length(tg_calc_result)
    tg_data_line = tg_calc_result{i}.data;
    if tg_calc_result{i}.is_ok
      if contains(tg_calc_result{i}.filename, 'ff') && length(tg_data_line) > 2
        tg_y = zeros(length(tg_data_line), 1);
        tg_x = tg_data_line{1}.result.x;
        tg_i = zeros(length(tg_x), length(tg_data_line));

        if contains(tg_calc_result{i}.filename, 'pp') || contains(tg_calc_result{i}.filename, 'fd')
          tg_fft_i = zeros(fix(length(tg_x) / 2 + 1), length(tg_data_line));
          if contains(tg_calc_result{i}.filename, 'fd')
              tg_i_1d = zeros(length(tg_data_line), 1);
          end
          for j = 1:1:length(tg_data_line)
            if contains(tg_calc_result{i}.filename, 'pp')
              tg_y(j) = tg_data_line{j}.sample.tg.pp;
            else
              tg_y(j) = 1 / tg_data_line{j}.probe.wfc;
              tg_i_1d(j) = sum(tg_data_line{j}.result.intensity.x(1:length(tg_x)));
            end
            tg_i(:, j) = tg_data_line{j}.result.intensity.x;
            [spectrum, width] = intensityResultSpectrum(tg_data_line{j}.result.intensity.x);
            tg_fft_i(:, j) = spectrum / max(spectrum);
          end
          if contains(tg_calc_result{i}.filename, 'fd')
            display1DImage(tg_i_1d, tg_y, 10);
            [tg_i_1d_fft, tg_i_1d_dw] = intensityResultSpectrum(tg_i_1d - mean(tg_i_1d));
            tg_i_1d_df = ((0:1:tg_i_1d_dw) / (max(tg_y) - min(tg_y) + tg_y(2) - tg_y(1)))';
            displayintensitySpectrum(tg_i_1d_fft, tg_i_1d_df, 11);
          end
          tg_f = ((0:1:width) / (max(tg_x) - min(tg_x) + tg_x(2) - tg_x(1)))';
          if contains(tg_calc_result{i}.filename, 'pp')
            y_label = 'Distance [mm]';
          else
            y_label = 'Curvature radius [mm]';
          end
          display2DImage(tg_fft_i, tg_f, tg_y, 'FFT image', 'f [mm^{-1}]', y_label, 12);
          saveas(gcf, append(tg_calc_result{i}.wd, 'ff-fft_2d_scale.png'));
          imwrite(tg_fft_i * 255, jet, append(tg_calc_result{i}.wd, 'ff-fft_2d.png'), 'png');
          tg_fft_i = convertScale(tg_fft_i, 1e10);
          display2DImage(tg_fft_i, tg_f, tg_y, 'FFT image', 'f [mm^{-1}]', y_label, 13);
          imwrite(tg_fft_i * 255, jet, append(tg_calc_result{i}.wd, 'ff-fft-log_2d.png'), 'png');
          saveas(gcf, append(tg_calc_result{i}.wd, 'ff-fft_log_2d_scale.png'));
        end

        if contains(tg_calc_result{i}.filename, 'ncb')
          for j = 1:1:length(tg_data_line)
            tg_y(j) = tg_data_line{j}.probe.intensity;
            tg_i(:, j) = tg_y(j) * tg_data_line{j}.result.intensity.x;
          end
          display2DImage(tg_i, tg_x, tg_y, 'FF image', 'x [\mu{m}]', 'I [a.u.]', 10);
          imwrite(tg_i * 255 / max(tg_y), jet, append(tg_calc_result{i}.wd, 'ff_2d.png'), 'png');
        end

        if contains(tg_calc_result{i}.filename, 'g2s')
          tg_i_1d = zeros(length(tg_data_line), 1);
          tg_fft_i = zeros(fix(length(tg_x) / 2 + 1), length(tg_data_line));
          for j = 1:1:length(tg_data_line)
            tg_y(j) = tg_data_line{j}.sys.g2s;
            tg_i(:, j) = tg_data_line{j}.result.intensity.x(1:length(tg_x));
            tg_i_1d(j) = sum(tg_data_line{j}.result.intensity.x(1:length(tg_x)));
            [spectrum, width] = intensityResultSpectrum(tg_data_line{j}.result.intensity.x);
            tg_fft_i(:, j) = spectrum / max(spectrum);
          end
          tg_f = ((0:1:width) / (max(tg_x) - min(tg_x) + tg_x(2) - tg_x(1)))';
          display1DImage(tg_i_1d, tg_y, 10);
          [tg_i_1d_fft, tg_i_1d_dw] = intensityResultSpectrum(tg_i_1d - mean(tg_i_1d));
          tg_i_1d_df = ((0:1:tg_i_1d_dw) / (max(tg_y) - min(tg_y) + tg_y(2) - tg_y(1)))';
          displayintensitySpectrum(tg_i_1d_fft, tg_i_1d_df, 11);
          calc_string = ["Z, mm", "I_x, a.u.", "f, mm^-1", "p, um", "FFT, a.u."];
          writematrix(calc_string, tg_calc_result{i}.filename, 'Sheet', 'TGSignalCalc', 'Range', 'A1');
          calc_string = [tg_y, tg_i_1d];
          writematrix(calc_string, tg_calc_result{i}.filename, 'Sheet', 'TGSignalCalc', 'Range', 'A2');
          calc_string = [tg_i_1d_df, (tg_i_1d_df.^-1) * 1000, tg_i_1d_fft];
          writematrix(calc_string, tg_calc_result{i}.filename, 'Sheet', 'TGSignalCalc', 'Range', 'C2');
          saveas(gcf, append(tg_calc_result{i}.wd, 'ff_fft_scale.png'));
          display2DImage(tg_fft_i, tg_f, tg_y, 'FFT image', 'f [mm^{-1}]', 'Distance [mm]', 12);
          saveas(gcf, append(tg_calc_result{i}.wd, 'ff-fft_2d_scale.png'));
          imwrite(tg_fft_i * 255, jet, append(tg_calc_result{i}.wd, 'ff-fft_2d.png'), 'png');
          tg_fft_i = convertScale(tg_fft_i, 1e10);
          display2DImage(tg_fft_i, tg_f, tg_y, 'FFT image', 'f [mm^{-1}]', 'Distance [mm]', 13);
          imwrite(tg_fft_i * 255, jet, append(tg_calc_result{i}.wd, 'ff-fft-log_2d.png'), 'png');
          saveas(gcf, append(tg_calc_result{i}.wd, 'ff-fft_log_2d_scale.png'));
        end

        if contains(tg_calc_result{i}.filename, 'ib')
          for j = 1:1:length(tg_data_line)
            tg_y(j) = tg_data_line{j}.sample.tg.beta;
            tg_i(:, j) = tg_data_line{j}.result.intensity.x(1:length(tg_x));
          end
        end

        if contains(tg_calc_result{i}.filename, 'g2psd') && ~contains(tg_calc_result{i}.filename, 'nosample')
          is_ref = false;
          for j = 1:1:length(tg_calc_result)
            if contains(tg_calc_result{j}.filename, 'g2psd') && contains(tg_calc_result{j}.filename, 'nosample') && ...
               tg_calc_result{j}.is_ok
              [tg_i_ref, ~, ~] = visualizeG2PSD(tg_calc_result{j}.data, tg_calc_result{j}.wd);
              is_ref = true;
            end
          end
          [tg_i, tg_y, tg_x] = visualizeG2PSD(tg_data_line, tg_calc_result{i}.wd);
          if is_ref
            tg_i_diff = tg_i_ref / max(max(tg_i_ref)) - tg_i / max(max(tg_i));
            tg_i_diff = tg_i_diff - min(min(tg_i_diff));
            tg_i_diff = tg_i_diff / max(max(tg_i_diff));
            display2DImage(tg_i_diff, tg_y, tg_x, 'FF difference image', 'x [{\mu}m]', 'Distance [mm]', 10);
            imwrite(tg_i_diff * 255, jet, append(tg_calc_result{i}.wd, 'ff_diff_2d.png'), 'png');
            tg_i_diff = convertScale(tg_i_diff, 1e1);
            imwrite(tg_i_diff * 255, jet, append(tg_calc_result{i}.wd, 'ff-diff-log_2d.png'), 'png');
            display2DImage(tg_i_diff, tg_y, tg_x, 'FF difference image', 'x [{\mu}m]', 'Distance [mm]', 11);
          end
        end

        if ~contains(tg_calc_result{i}.filename, 'g2psd')
          display2DImage(tg_i, tg_x, tg_y, 'FF image', 'x [mm]', 'Distance [mm]', 8);
          saveas(gcf, append(tg_calc_result{i}.wd, 'ff_2d_scale.png'));
          imwrite(tg_i * 255, jet, append(tg_calc_result{i}.wd, 'ff_2d.png'), 'png');
          tg_i = convertScale(tg_i, 1e5);
          imwrite(tg_i * 255, jet, append(tg_calc_result{i}.wd, 'ff-log_2d.png'), 'png');
          display2DImage(tg_i, tg_x, tg_y, 'FF image', 'x [mm]', 'Distance [mm]', 9);
          saveas(gcf, append(tg_calc_result{i}.wd, 'ff_log_2d_scale.png'));
        end
      end
    end
  end
end

function displayintensitySpectrum(spectrum, f, i)
  figure(i); plot((f.^-1) * 1000, spectrum);
  xlabel('Spatial period, {\mu}m');
  ylabel('S, arb.un.');
  yscale log; grid on;
  ylim([1E-5 1]);
  title('Intensity Spectrum');
end

function [spectrum, width] = intensityResultSpectrum(intensity)
  cn = fft(intensity);
  spectrum = abs(cn).^2;
  spectrum = spectrum / max(spectrum);
  width = fix(length(spectrum) / 2);
  spectrum = spectrum(1: width + 1);
end

function [tg_i, tg_g2psd, tg_x] = visualizeG2PSD(tg_data_line, wd)
  tg_g2psd = zeros(length(tg_data_line), 1);
  tg_x = tg_data_line{1}.result.x;
  tg_i = zeros(length(tg_x), length(tg_data_line));
  for j = 1:1:length(tg_data_line)
    tg_g2psd(j) = tg_data_line{j}.sys.ff;
    tg_i(:, j) = tg_data_line{j}.result.intensity.x;
  end
  display2DImage(tg_i, tg_x, tg_g2psd, 'FF image', 'x [{\mu}m]', 'Distance [mm]', 8);
  imwrite(tg_i * 255, jet, append(wd, 'ff_2d.png'), 'png');
  tg_i_log = convertScale(tg_i, 1);
  imwrite(tg_i_log * 255, jet, append(wd, 'ff-log_2d.png'), 'png');
  display2DImage(tg_i_log, tg_x, tg_g2psd, 'FF image', 'x [{\mu}m]', 'Distance [mm]', 9);
end

function display2DImage(image, x, y, str, x_label, y_label, i)
  figure(i); imagesc(y, x, image); colormap(jet); colorbar; title(str); ylabel(x_label); xlabel(y_label);
end

function display1DImage(intensity, x, i)
  figure(i); plot(x, intensity); xlabel('x, mm');
  ylabel('I, arb.un.');
  grid on;
  xlim([min(x) max(x)]);
  ylim([min(intensity) max(intensity)]);
  title('Intensity');
end

function display1DError(error, n, ymin, i)
  error = error / max(error);
  figure(i);
  plot(n, error, '-o', 'MarkerSize', 8, 'MarkerEdgeColor','blue', ...
       'MarkerFaceColor', [.5 .5 1], 'Color', 'blue');
  xlabel('Line#');
  ylabel('Error, a.u.');
  yscale log;
  grid on;
  xlim([min(n) max(n)]);
  ylim([ymin 1]);
end

function image = convertScale(image, image_coef)
  image = log10(image_coef * image + 1);
  image = image / max(max(image));
end