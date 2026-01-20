function [is_ok, x, intensity, message, io, pump, probe, sample] = processTGCalculation(grating, pump, probe, sample, sys, io)
  x = 0;
  intensity.x = 0;
  intensity.y = 0;
  is_ok = false;
  if isnan(grating.p) || isnan(pump.lambda) || isnan(probe.lambda) || isnan(sys.g2s) || ...
     isnan(sys.ff) || isnan(sys.psd) || isnan(sys.ds)
    message = 'Wrong data: missing mandatory value for grating period, pump or probe wavelenght, probe intensity, distances, image aperture or step';
    warning(append('  Line #', io.i, '. ', message), 'Data Range warning');
  elseif grating.p <= 0 || pump.waist < 0 || pump.lambda <= 0 || pump.intensity < 0 || ...
         probe.waist < 0 || probe.lambda <= 0 || probe.intensity < 0 || sample.t < 0 || ...
         sys.g2s <= 0 || sys.ff <= 0 || sys.psd <= 0 || sys.ds <= 0
    message = 'Wrong data: only positive values are allowed for grating period; beams intensities, waists and wavelengths; distances, image aperture and step';
    warning(append('  Line #', io.i, '. ', message), 'Data Range warning');
  elseif (~isnan(pump.wfc) && 2 * grating.p^2 / pump.lambda * abs(pump.wfc) > 1) || ...
         (~isnan(probe.wfc) && 2 * grating.p^2 / probe.lambda * abs(probe.wfc) > 1)
    message = 'Wrong data: approximation is not possible since wavefront curvature is too big';
    warning(append('  Line #', io.i, '. ', message), 'Data Range warning');
  elseif isnan(probe.intensity) || probe.intensity == 0
    message = 'Wrong data: probe beam intensity is not set';
    warning(append('  Line #', io.i, '. ', message), 'Data Range warning');
  elseif sys.psd < sys.ds
    message = 'Wrong data: PSD aperture is less then calculation step';
    warning(append('  Line #', io.i, '. ', message), 'Data Range warning');
  elseif grating.df >= 1 || grating.df < 0
    message = 'Wrong data: grating duty factor must belong to the range [0, 1)';
    warning(append('  Line #', io.i, '. ', message), 'Data Range warning');
  elseif sys.dsdp < 1
    message = 'Wrong data: division parameter is less than 1';
    warning(append('  Line #', io.i, '. ', message), 'Data Range warning');
  elseif ~contains(grating.slit, 'cos') && ~contains(grating.slit, 'square') && ~contains(grating.slit, 'delta')
    message = 'Wrong configuration: three grating types are considred: cos, delta or square';
    warning(append('  Line #', io.i, '. ', message), 'Configuration warning');
  else
    if isnan(grating.df)
      grating.df = 0;
    end
    if isnan(pump.pd)
      pump.pd = 0;
    end
    if isnan(pump.angle)
      pump.angle = 0;
    end
    if isnan(pump.wfc)
      pump.wfc = 0;
    end
    if isnan(pump.waist)
      pump.waist = 0;
    end
    if isnan(pump.intensity)
      pump.intensity = 0;
    end
    if isnan(probe.pd)
      probe.pd = 0;
    end
    if isnan(probe.angle)
      probe.angle = 0;
    end
    if isnan(probe.wfc)
      probe.wfc = 0;
    end
    if isnan(probe.waist)
      probe.waist = 0;
    end
    if isnan(sample.t)
      sample.t = 0;
    end
    % it is considered no absorbance for now -> imag(sample.n2) = 0
    sample.n2 = complex(real(sample.n2), 0);
    if isnan(real(sample.n2))
      sample.n2 = complex(0, 0);
    end
    if isnan(sys.dsdp) || isinf(sys.dsdp)
      sys.dsdp = fix(sys.ds / grating.p) + 1;
    end
    [is_ok, x, intensity, message, pump, probe, sample, io] = calculateTGFFImage(grating, pump, probe, sample, sys, io);
  end
end

function [is_ok, x, intensity, message, pump, probe, sample, io] = calculateTGFFImage(grating, pump, probe, sample, sys, io)
  if strlength(io.filename) == 0
    io.sdp = false;
  end
  if io.sdp
    io.line = append('LineN', io.i);
  end
  x = 0;
  intensity.x = 0;
  intensity.y = 0;
  sys.points = fix(sys.psd / sys.ds) + 1;
  if mod(sys.points, 2) == 0
    sys.points = sys.points + 1;
  end
  sys.points_nodsdp = sys.points;
  sys.psd = sys.ds * sys.points;
  sys.points = sys.points * sys.dsdp;
  edges = 10.^linspace(-20, 0, 21);

  [pump.cn, probe.cn] = calculateCn(grating, pump, probe, edges, io);

  [is_ok, sample, message] = calculateAm(grating, pump, probe, sys, edges, sample, io);

  if is_ok
    [x, intensity] = calculateOutpuIntensity(probe, grating, sample, sys, io);
  end
end

function [x, intensity] = calculateOutpuIntensity(probe, grating, sample, sys, io)
  [B, D] = calculateBD(probe.lambda, grating.p, probe.angle);
  x = linspace(-sys.psd / 2, sys.psd / 2, sys.points)';
  [A, C, E, F, G] = calculateACEFG(probe.waist, probe.lambda, sys.ff, ...
                                   probe.wfc, sys.g2s, sample.aperture);
  [~, intensity.x] = diffractionGS1DFinite(probe.cn, sample.am, x, A, B, C, D, E, F, G);
  if probe.waist ~= 0
    intensity.y = diffractionG01DFinite(sys.ff, probe.wfc, probe.waist, probe.lambda, x);
  else
    intensity.y = ones(1, length(x));
  end
  if sys.dsdp == 1
    x_out = x;
    i_out_x = intensity.x;
    i_out_y = intensity.y;
  else
    x_out = linspace(-sys.psd / 2, sys.psd / 2, sys.points_nodsdp)';
    i_out_x = zeros(sys.points_nodsdp, 1);
    i_out_y = zeros(sys.points_nodsdp, 1);
    for i = 1:1:sys.points_nodsdp
      i_out_x(i) = sum(intensity.x((i - 1) * sys.dsdp + 1:i * sys.dsdp)) / sys.dsdp;
      i_out_y(i) = sum(intensity.y((i - 1) * sys.dsdp + 1:i * sys.dsdp)) / sys.dsdp;
    end
  end
  intensity.x = i_out_x;
  intensity.y = i_out_y;
  x = x_out;

  display1DImage(intensity.x, x, 3);
  if probe.waist ~= 0
    displayTwo1DImages(intensity, x, x', 3);
  end
  if io.sdp
    saveas(gcf, append(io.wd, 'i_', io.line, '.png'));
  end

  if length(x) * length(x) < 1E10
    diffraction_image = intensity.x' .* intensity.y;
    diffraction_image_visual = diffraction_image / max(max(diffraction_image));
    display2DImage(diffraction_image_visual, x, x', 'Intensity of Diffracted Wave', 4);
    if io.ddp
      [i_x_fft, i_x_f] = intensityResultSpectrum(intensity.x - mean(intensity.x), x);
      displayIntensityLinSpectrum(i_x_fft, i_x_f, 17);
      if io.sdp
        saveas(gcf, append(io.wd, 'i_x_fft_scale_', io.line, '.png'));
      end
      displayIntensitySpectrum(i_x_fft, i_x_f, 18);
      if io.sdp
        saveas(gcf, append(io.wd, 'i_x_fft_log_scale_', io.line, '.png'));
      end
    end
    if io.sdp
      imwrite(diffraction_image * 255, jet, append(io.wd, 'i_2d_', io.line, '.png'), 'png');
      saveas(gcf, append(io.wd, 'i_2d_scale_', io.line, '.png'));
    end
    [diffraction_image_visual, image_coef] = convertScale(diffraction_image_visual);
    display2DImage(diffraction_image_visual, x, x', append('Intensity of Diffracted Wave LOG(', string(image_coef), ')'), 5);
    if io.sdp
      imwrite(diffraction_image_visual * 255, jet, append(io.wd, 'i_2d_log_', io.line, '.png'), 'png');
      saveas(gcf, append(io.wd, 'i_2d_log_scale_', io.line, '.png'));
    end
  end
  if strlength(io.filename) ~= 0 && io.sdp
    calc_string = ["x, mm", "I_x, a.u.", "y, mm", "I_y, a.u."];
    writematrix(calc_string, io.filename, 'Sheet', io.line, 'Range', 'K1');
    calc_string = [x, intensity.x, x, intensity.y];
    writematrix(calc_string, io.filename, 'Sheet', io.line, 'Range', 'K2');
  end
end

function [image, image_coef] = convertScale(image)
  image_coef = 1 / min(min(image));
  if image_coef > 1E5
    image_coef = 1E5;
  end
  image = log10(image_coef * image + 1);
  image = image / max(max(image));
end

function intensity = diffractionG01DFinite(ff, wfc, waist, lambda, y)
  ff = ff / (1 + ff * wfc);
  A = -2 * pi^2 * (ff * wfc - 1)^2 / (lambda * ff * complex(lambda * ff * waist^2, -2 * pi));
  u = exp(A * y.^2);
  intensity = abs(u).^2;
end

function display1DImage(intensity, x, i)
  intensity = intensity / max(intensity);
  figure(i); plot(x, intensity); xlabel('x, mm');
  ylabel('I, arb.un.');
  yscale log;
  grid on;
  xlim([min(x) max(x)]);
  ylim([1E-5 1]);
  title('Intensity');
end

function displayTwo1DImages(intensity, x, y, i)
  intensity.x = intensity.x / max(intensity.x);
  intensity.y = intensity.y / max(intensity.y);
  figure(i);
  plot(x, intensity.x);
  hold on;
  plot(y, intensity.y);
  hold off;
  xlabel('x, mm');
  ylabel('I, arb.un.');
  yscale log;
  grid on;
  ylim([1E-5 1]);
  xlim([min(x) max(x)]);
  title('Intensity');
  legend('0X','0Y');
end

function [is_ok, sample, message] = calculateAm(grating, pump, probe, sys, edges, sample, io)
  is_ok = true;
  message = '';
  beta = 2 * pi * sample.t * (probe.intensity + pump.intensity) * sample.n2 / probe.lambda;
  am = zeros(1, 2);
  am(2) = 1;
  pp = pump.intensity / (probe.intensity + pump.intensity);
  if pp < 0.01
    if probe.waist == 0 && probe.wfc == 0
      aperture = grating.p;
    elseif probe.waist == 0 && probe.wfc ~= 0
      aperture = probe.lambda * grating.p * (1 + sys.g2s * probe.wfc) / ...
                 (probe.lambda * (1 + sys.g2s * probe.wfc) - ...
                 2 * grating.p^2 * probe.wfc);
    else
      aperture = 4 * (1 + sys.g2s * abs(probe.wfc)) / probe.waist;
    end
  elseif pp > 0.99
    if pump.waist == 0 && pump.wfc == 0
      aperture = grating.p;
    elseif pump.waist == 0 && pump.wfc ~= 0
      aperture = pump.lambda * grating.p * (1 + sys.g2s * pump.wfc) / ...
                 (pump.lambda * (1 + sys.g2s * pump.wfc) - ...
                 2 * grating.p^2 * pump.wfc);
    else
      aperture = 4 * (1 + sys.g2s * abs(probe.wfc)) / probe.waist;
    end
  else
    if probe.waist == 0 
      aperture = 1.5 * sys.psd;
    else
      aperture = 4 * (1 + sys.g2s * abs(probe.wfc)) / probe.waist;
    end
  end
  if beta ~= 0
    x = -aperture / 2:grating.p / 50:aperture / 2;
    if mod(length(x), 2) == 0
      x = linspace(-aperture / 2, aperture / 2, length(x) - 1);
    end
    x = x';
    intensity.pump = calculateSIntensity(grating.p, pump.lambda, sys.g2s, ...
                                              pump.wfc, pump.waist, pump.angle, pump.cn, x);
    intensity.pump = pp * intensity.pump / max(intensity.pump);
    intensity.probe = calculateSIntensity(grating.p, probe.lambda, sys.g2s, ...
                                                     probe.wfc, probe.waist, probe.angle, probe.cn, x);
    intensity.probe = (1 - pp) * intensity.probe / max(intensity.probe);
    intensity.sum = intensity.pump + intensity.probe;
    if io.ddp
      display1DIntensities(intensity, x);
      if io.sdp
        saveas(gcf, append(io.wd, 'i_beta_', io.line, '.png'));
      end
    end
    point_ns = 0;
    for i = 1:1:length(x)
      if x(i) >= -grating.p * 5 && point_ns == 0
        point_ns = i;
      end
      if x(i) > grating.p * 5
        break;
      end
    end
    sample.tg.x = x(point_ns:i - 1);
    sample.tg.pd = intensity.sum(point_ns:i - 1);
    sample.tg.beta = beta;
    sample.tg.pp = pp;
    if strlength(io.filename) ~= 0 && io.sdp
      calc_string = ["beta, rad/I", "x, mm", "I_beta, a.u."];
      writematrix(calc_string, io.filename, 'Sheet', io.line, 'Range', 'G1');
      writematrix(sample.tg.beta, io.filename, 'Sheet', io.line, 'Range', 'G2');
      calc_string = [sample.tg.x, sample.tg.pd];
      writematrix(calc_string, io.filename, 'Sheet', io.line, 'Range', 'H2');
    end
    intensity.sum = intensity.sum * beta;
    tr = complex(cos(intensity.sum), sin(intensity.sum));
    am = fftshift(fft(tr));
    spectrum = abs(am).^2;
    am = am / max(sqrt(spectrum));
    spectrum = spectrum / max(spectrum);
    width = fix(length(spectrum) / 2);
    am = fourierCoefficientUpdate(am, width, spectrum, edges);
    size_am = size(am);
    if (size_am(1) == 1 && size_am(2) == 2) || (size_am(1) == 2 && size_am(2) == 1)
      message = 'Wrong input data: sample phase effectivity is set as zero. Please check input parameters.';
      is_ok = false;
    end
  end
  sample.aperture = aperture;
  sample.am = am;
end

function intensity = calculateSIntensity(p, lambda, g2s, wfc, waist, angle, cn, x)
  [B, D] = calculateBD(lambda, p, angle);
  [A, C] = calculateAC(waist, lambda, g2s, wfc);
  [~, intensity] = diffraction1DFinite(cn, A, B, C, D, x);
end

function [u, intensity] = diffraction1DFinite(cn, A, B, C, D, x)
  n = D * real(cn(:, 1));
  x = C * x + B;
  u = zeros(length(x), 1);
  for i = 1:1:length(x)
    for j = 1:1:length(n)
      exp_a = exp(-A * (x(i) + n(j))^2);
      u(i) = u(i) + exp_a * cn(j, 2);
    end
  end
  intensity = abs(u).^2;
end

function [u, intensity] = diffractionGS1DFinite(cn, am, x, A, B, C, D, E, F, G)
  u = zeros(length(x), 1);
  xG = G * x;
  x = C * x + B;
  n = D * real(cn(:, 1));
  m = E * real(am(:, 1));
  for i = 1:1:length(x)
    for j = 1:1:length(n)
      for k = 1:1:length(m)
        c = cn(j, 2) * am(k, 2);
        exp_a = exp(F * (m(k) / 2 - xG(i))^2 - A * (x(i) + n(j) + m(k))^2);
        u(i) = u(i) + exp_a * c;
      end
    end
  end
  intensity = abs(u).^2;
  if max(intensity) ~= 0
    u = u / sqrt(max(intensity));
  end
end

function [B, D] = calculateBD(lambda, p, angle)
  B = - sin(angle) / (2 * lambda);
  D = 1 / p;
end

function [A, C] = calculateAC(waist, lambda, distance, wfc)
  A = pi^2 * 2 * lambda * distance / complex(lambda * distance * waist^2, -2 * pi * (1 + wfc * distance));
  C = -1 / (lambda * distance);
end

function [A, C, E, F, G] = calculateACEFG(waist, lambda, ff, wfc, g2s, aperture)
  L = ff - g2s;
  A = pi^2 * 2 * lambda * ff / complex(lambda * ff * waist^2, -2 * pi * (1 + wfc * ff));
  C = -1 / (lambda * ff);
  E = L / (aperture * ff);
  F = pi * lambda * L * g2s / ff * complex(0, 1);
  G = 1 / (lambda * ff);
end

function cn_u = fourierCoefficientUpdate(cn, width, spectrum, edges)
  [N, ~] = histcounts(spectrum, edges);
  index = 1;
  for i = length(N) - 1:-1:2
    if index == 1 && sum(N(i:length(N))) > 40
      index = i;
    end
  end
  cn_u = zeros(length(cn), 2);
  i = 1;
  for f = -width:1:width
    n = f + width + 1;
    if spectrum(n) > edges(index)
      cn_u(i, 1) = f;
      cn_u(i, 2) = cn(n);
      i = i + 1;
    end
  end
  cn_u = cn_u(1:i - 1, :);
end

function displayTwoSpectra(width_pump, spectrum_pump, width_probe, spectrum_probe)
  f0 = -width_pump:1:width_pump;
  f1 = -width_probe:1:width_probe;
  figure(1);
  bar(f0, spectrum_pump);
  hold on;
  bar(f1, spectrum_probe);
  hold off;
  xlabel('n');
  ylabel('S, arb.un.');
  title('Grating Cn');
  legend('Beam pump','Beam probe');
end

function display1DIntensities(intensity, x)
  figure(2);
  plot(x, intensity.pump);
  hold on;
  plot(x, intensity.probe);
  plot(x, intensity.sum);
  hold off;
  xlabel('x, mm');
  ylabel('I, arb.un.');
  yscale log;
  grid on;
  xlim([min(x) max(x)]);
  title('Intensities on the sample');
  legend('Beam pump','Beam probe', 'Sum');
end

function [pump_cn, probe_cn] = calculateCn(grating, pump, probe, edges, io)
  if ~contains(grating.slit, 'phase')
    [spectrum, cn_re, width] = calculateCnRe(grating);
    cn_im = zeros(length(cn_re), 1);
    pump_cn = complex(cn_re, cn_im);
    probe_cn = complex(cn_re, cn_im);
    width_pump = width;
    width_probe = width;
    spectrum_pump = spectrum;
    spectrum_probe = spectrum;
  else
    [spectrum_pump, pump_cn, width_pump] = calculateCnComplex(grating, pump.lambda, pump.pd);
    [spectrum_probe, probe_cn, width_probe] = calculateCnComplex(grating, probe.lambda, probe.pd);
  end
  if io.ddp
    displayTwoSpectra(width_pump, spectrum_pump, width_probe, spectrum_probe);
    if io.sdp
      saveas(gcf, append(io.wd, 'cn_coef_', io.line, '.png'));
    end
  end
  pump_cn = fourierCoefficientUpdate(pump_cn, width_pump, spectrum_pump, edges);
  probe_cn = fourierCoefficientUpdate(probe_cn, width_probe, spectrum_probe, edges);
  if strlength(io.filename) ~= 0 && io.sdp
    calc_string = ["n_pump", "Re(pump.cn)", "Im(pump.cn)", "n_probe", "Re(probe.cn)", "Im(probe.cn)"];
    writematrix(calc_string, io.filename, 'Sheet', io.line, 'Range', 'A1');
    calc_string = [real(pump_cn(:, 1)), real(pump_cn(:, 2)), imag(pump_cn(:, 2))];
    writematrix(calc_string, io.filename, 'Sheet', io.line, 'Range', 'A2');
    calc_string = [real(probe_cn(:, 1)), real(probe_cn(:, 2)), imag(probe_cn(:, 2))];
    writematrix(calc_string, io.filename, 'Sheet', io.line, 'Range', 'D2');
  end
end

function [spectrum, cn_re, width] = calculateCnRe(grating)
  if contains(grating.slit, 'delta') || (contains(grating.slit, 'square') && grating.df == 0)
    width = 300;
    cn_re = ones(2 * width + 1, 1) / (2 * width + 1);
  elseif contains(grating.slit, 'cos')
    width = 1;
    cn_re = 0.25 * ones(2 * width + 1, 1);
    cn_re(width + 1) = 0.5;
  else
    % Square grating.slit type is considered
    width0 = 300;
    cn_re = zeros(2 * width0 + 1, 1);
    c0 = grating.df;
    c1 = pi * grating.df;
    cn_re(width0 + 1) = c0;
    for i = 1:1:width0
      x = c1 * i;
      cn_re(width0 + 1 + i) = c0 * sin(x) / x;
      cn_re(width0 + 1 - i) = cn_re(width0 + 1 + i);
    end
    spectrum = cn_re.^2 / c0^2;
    i = 1;
    while spectrum(i) < 1E-4
      i = i + 1;
    end
    i = i - 1;
    width = width0 - i;
    cn_re = cn_re(width0 - width + 1:width0 + width + 1);
  end
  spectrum = cn_re.^2;
end

function [spectrum, cn, width] = calculateCnComplex(grating, lambda, pd)
  cos_phi = cos(pd);
  sin_phi = sin(pd);
  if contains(grating.slit, 'delta') || (contains(grating.slit, 'square') && grating.df == 0)
    width0 = 300;
    pit_re = zeros(2 * width0 + 1, 1);
    pit_im = zeros(2 * width0 + 1, 1);
    pit_re(1) = 1;
    pit_re(2 * width0 + 1) = 1;
    for i = 2:1:2 * width0
      pit_re(i) = cos_phi;
      pit_im(i) = sin_phi;
    end
  else
    width0 = fix(grating.p / lambda) + 1;
    pit_re = ones(2 * width0 + 1, 1);
    pit_im = zeros(2 * width0 + 1, 1);
    if contains(grating.slit, 'cos')
      x = linspace(0, grating.p, 2 * width0 + 2);
      for i = 1:1:2 * width0 + 1
        phase = pd * cos(2 * pi / grating.p * x(i));
        pit_re(i) = cos(phase);
        pit_im(i) = sin(phase);
      end
    else
      points = fix(grating.p * grating.df / lambda);
      for i = 0:1:points
        pit_re(1 + i) = cos_phi;
        pit_re(2 * width0 + 1 - i) = cos_phi;
        pit_im(1 + i) = sin_phi;
        pit_im(2 * width0 + 1 - i) = sin_phi;
      end
    end
  end
  pit = complex(pit_re, pit_im);
  cn = fftshift(fft(pit));
  spectrum = abs(cn).^2;
  cn = cn / max(sqrt(spectrum));
  spectrum = spectrum / max(spectrum);
  i = 1;
  while spectrum(i) < 1E-4
    i = i + 1;
  end
  i = i - 1;
  width = width0 - i;
  cn = cn(width0 - width + 1:width0 + width + 1);
  spectrum = spectrum(width0 - width + 1:width0 + width + 1);
end

function display2DImage(image, x, y, str, i)
  figure(i); imagesc(y, x, image); colormap(jet); colorbar; title(str);
end

function [spectrum, f] = intensityResultSpectrum(intensity, x)
  cn = fft(intensity);
  spectrum = abs(cn).^2;
  spectrum = spectrum / max(spectrum);
  width = fix(length(spectrum) / 2);
  spectrum = spectrum(1: width + 1);
  f = ((0:1:width) / (max(x) - min(x) + x(2) - x(1)));
end

function displayIntensitySpectrum(spectrum, f, i)
  figure(i); 
  plot(f, spectrum);
  xlabel('Spatial frequency, mm^{-1}');
  ylabel('S, arb.un.');
  yscale log; grid on;
  ylim([1E-5 1]);
  xlim([min(f) max(f)]);
  title('Intensity Spectrum');
end

function displayIntensityLinSpectrum(spectrum, f, i)
  figure(i); plot(f, spectrum);
  xlabel('Spatial frequency, {\mu}m');
  ylabel('S, arb.un.');
  grid on;
  ylim([0 1]);
  xlim([min(f) max(f)]);
  title('Intensity Spectrum');
end
