%%Function For Shape Generation Via the Supefomula

function SuperFormula3D(a1, b1, m1, n11, n21, n31, a2, b2, m2, n12, n22, n32)
    % Benchmark if not values provides
    if nargin < 12
        a1 = 1; b1 = 1; m1 = 6; n11 = 1; n21 = 1; n31 = 1; %Keep a's and b's equal to 1 for not weird shapes 
        a2 = 1; b2 = 1; m2 = 6; n12 = 1; n22 = 1; n32 = 1; %Kepp a's and b's equal to 1 for not weirid shapes 
    end

    % Other Parameters (don't change)
    rlong = 1; rlat = 1;

    % Create grid of phi and theta values (increase value for higher mesh -
    % Not neccesary)
    steps = 100; %Panta 100
    [phi, theta] = meshgrid(linspace(-pi/2, pi/2, steps), linspace(-pi, pi, steps)); %Jin's Parameters

    % Evaluation of r1 and r2 based on inputs 
    r1 = ssf(theta, a1, b1, m1, n11, n21, n31);
    r2 = ssf(phi, a2, b2, m2, n12, n22, n32);

    % Evaluation of x, y, and z based on Inputs 
    x = rlong * r1 .* cos(theta) .* r2 .* cos(phi);
    y = rlong * r1 .* sin(theta) .* r2 .* cos(phi);
    z = rlat * r2 .* sin(phi);

    % For Figure
    figure('Color', 'w', 'Position', [100, 100, 800, 600]);
    surf(x, y, z, 'EdgeColor', 'none', 'FaceColor', 'interp', 'FaceLighting', 'gouraud');
    light('Position', [1 1 1], 'Style', 'infinite');
    lighting gouraud
    view(45, 30);
    axis equal tight off;
    title_str = sprintf('3D Supershape\n');
    title_str = [title_str sprintf('a1=%.2f, b1=%.2f, m1=%.2f, n11=%.2f, n21=%.2f, n31=%.2f\n', a1, b1, m1, n11, n21, n31)];
    title_str = [title_str sprintf('a2=%.2f, b2=%.2f, m2=%.2f, n12=%.2f, n22=%.2f, n32=%.2f', a2, b2, m2, n12, n22, n32)];
    title(title_str, 'FontSize', 10);

    % Massage for accepting or declining a geometry   
     user_response = input('Do you want to save this shape? (y/n): ', 's');
    if strcmpi(user_response, 'y')
        % Apothikefse sto directoy 
        save_dir = 'C:\Users\andre\OneDrive\Υπολογιστής\ImperialCollage\MATE70004-ResearchProject\Shapes\FiguresForReport';
        if ~exist(save_dir, 'dir')
            mkdir(save_dir);
        end
        existing_files = dir(fullfile(save_dir, 'Shape*.png'));

        % Number of shape and 
        num_existing_shapes = numel(existing_files) + 70;
        new_filename = sprintf('Shape%d.png', num_existing_shapes + 1);
        full_path_png = fullfile(save_dir, new_filename);
        saveas(gcf, full_path_png);

        % Create JSON file and save it in the working directoy
        json_filename = sprintf('input_shape%d.json', num_existing_shapes + 1);
        full_path_json = fullfile(save_dir, json_filename);
        %Jin's Parameters for the input file
        json_str = sprintf('{%s', newline);
        json_str = [json_str, '    "shape":"Superformula3D",', newline];
        json_str = [json_str, sprintf('    "shape_kwargs":{"rho":10.0,"m":%.1f,"n1":%.1f,"n2":%.1f,"n3":%.1f,"a":%.1f,"b":%.1f,"m2":%.1f,"n4":%.1f,"n5":%.1f,"n6":%.1f,"a2":%.1f,"b2":%.1f},', m1, n11, n21, n31, a1, b1, m2, n12, n22, n32, a2, b2), newline];
        json_str = [json_str, '    "materials_kwargs":{"cells":["Au 4.150925105415038"],', newline];
        json_str = [json_str, '    "SK_dict":{"Au":[-0.8759354,2.1869271,-0.2947458,-0.6343078,0.2991852,-0.0313741,1.2916726,-0.6334431,-0.8737013,0.2205823,0.0119527,0.0084755,0.0014077,0.0011205,-0.0025413,0.0011114,-0.0094128,0.0121376,-0.0122429,-0.0087764,8.9254145,17.7664004,4.6939335,4.1509251]}', newline];
        json_str = [json_str, '    },', newline];
        json_str = [json_str, '    "calculation":"HCG","construct_ase_nanoparticle_kwargs":{"offset":[[1e-6,1e-6,1e-6]]},', newline];
        json_str = [json_str, '    "calculator":"KPM","phi_calculator":"comsol","Phi_file":"Phi.npz",', newline];
        json_str = [json_str, '    "A":11.525177810989502,"B":14,"NNR":1500,"Nk":5000,"omega":2.4,"seednumber":114514,', newline];
        json_str = [json_str, '    "Ef":8.2226,"T":300,"sigma":0.06,"gamma":0.05,', newline];
        json_str = [json_str, '    "Nbatch":2,', newline];
        json_str = [json_str, '    "max_time":252000', newline];
        json_str = [json_str, '}', newline];
        fid = fopen(full_path_json, 'w');
        if fid == -1, error('Cannot create JSON file'); end
        fwrite(fid, json_str, 'char');
        fclose(fid, json_str, 'char');
        disp(['The 3D Supershape has been plotted and saved as "', full_path_png, '"']);%Ta shapes ginounte overwrite
        disp(['The corresponding input file has been saved as "', full_path_json, '"']);
    else
        disp('The shape was not saved.');%
    end
end

function r = ssf(theta, a, b, m, n1, n2, n3)
    t1 = abs(cos(m * theta / 4) / a).^n2;
    t2 = abs(sin(m * theta / 4) / b).^n3;
    r = (t1 + t2).^(-1/n1);
end