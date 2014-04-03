% Copyright (c) 2010, Oliver Woodford
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.
%EXPORT_FIG  Exports figures suitable for publication
%
% Examples:
%   im = export_fig
%   [im alpha] = export_fig
%   export_fig filename
%   export_fig filename -format1 -format2
%   export_fig ... -nocrop
%   export_fig ... -native
%   export_fig ... -m<val>
%   export_fig ... -r<val>
%   export_fig ... -a<val>
%   export_fig ... -q<val>
%   export_fig ... -<renderer>
%   export_fig ... -<colorspace>
%   export_fig ... -append
%   export_fig(..., handle)
%
% This function saves a figure or single axes to one or more vector and/or
% bitmap file formats, and/or outputs a rasterized version to the
% workspace, with the following properties:
%   - Figure/axes reproduced as it appears on screen
%   - Cropped borders (optional)
%   - Embedded fonts (vector formats)
%   - Improved line and grid line styles
%   - Anti-aliased graphics (bitmap formats)
%   - Render images at native resolution (optional for bitmap formats)
%   - Transparent background supported (pdf, eps, png)
%   - Semi-transparent patch objects supported (png only)
%   - RGB, CMYK or grayscale output (CMYK only with pdf, eps, tiff)
%   - Variable image compression, including lossless (pdf, eps, jpg)
%   - Optionally append to file (pdf, tiff)
%   - Vector formats: pdf, eps
%   - Bitmap formats: png, tiff, jpg, bmp, export to workspace 
%   
% This function is especially suited to exporting figures for use in
% publications and presentations, because of the high quality and
% portability of media produced.
%
% Note that the background color and figure dimensions are reproduced
% (the latter approximately, and ignoring cropping & magnification) in the
% output file. For transparent background (and semi-transparent patch
% objects), set the figure (and axes) 'Color' property to 'none'; pdf, eps
% and png are the only file formats to support a transparent background,
% whilst the png format alone supports transparency of patch objects.
%
% The choice of renderer (opengl, zbuffer or painters) has a large impact
% on the quality of output. Whilst the default value (opengl for bitmaps,
% painters for vector formats) generally gives good results, if you aren't
% satisfied then try another renderer.  Notes: 1) For vector formats (eps,
% pdf), only painters generates vector graphics. 2) For bitmaps, only
% opengl can render transparent patch objects correctly. 3) For bitmaps,
% only painters will correctly scale line dash and dot lengths when
% magnifying or anti-aliasing. 4) Fonts may be substitued with Courier when
% using painters.
%
% When exporting to vector format (pdf & eps) and bitmap format using the
% painters renderer, this function requires that ghostscript is installed
% on your system. You can download this from:
%   http://www.ghostscript.com
% When exporting to eps it additionally requires pdftops, from the Xpdf
% suite of functions. You can download this from:
%   http://www.foolabs.com/xpdf
%
%IN:
%   filename - string containing the name (optionally including full or
%              relative path) of the file the figure is to be saved as. If
%              a path is not specified, the figure is saved in the current
%              directory. If no name and no output arguments are specified,
%              the default name, 'export_fig_out', is used. If neither a
%              file extension nor a format are specified, a ".png" is added
%              and the figure saved in that format.
%   -format1, -format2, etc. - strings containing the extensions of the
%                              file formats the figure is to be saved as.
%                              Valid options are: '-pdf', '-eps', '-png',
%                              '-tif', '-jpg' and '-bmp'. All combinations
%                              of formats are valid.
%   -nocrop - option indicating that the borders of the output are not to
%             be cropped.
%   -m<val> - option where val indicates the factor to magnify the
%             on-screen figure dimensions by when generating bitmap
%             outputs. Default: '-m1'.
%   -r<val> - option val indicates the resolution (in pixels per inch) to
%             export bitmap outputs at, keeping the dimensions of the
%             on-screen figure. Default: sprintf('-r%g', get(0,
%             'ScreenPixelsPerInch')). Note that the -m and -r options
%             change the same property.
%   -native - option indicating that the output resolution (when outputting
%             a bitmap format) should be such that the vertical resolution
%             of the first suitable image found in the figure is at the
%             native resolution of that image. To specify a particular
%             image to use, give it the tag 'export_fig_native'. Notes:
%             This overrides any value set with the -m and -r options. It
%             also assumes that the image is displayed front-to-parallel
%             with the screen. The output resolution is approximate and
%             should not be relied upon. Anti-aliasing can have adverse
%             effects on image quality (disable with the -a1 option).
%   -a1, -a2, -a3, -a4 - option indicating the amount of anti-aliasing to
%                        use for bitmap outputs. '-a1' means no anti-
%                        aliasing; '-a4' is the maximum amount (default).
%   -<renderer> - option to force a particular renderer (painters, opengl
%                 or zbuffer) to be used over the default: opengl for
%                 bitmaps; painters for vector formats.
%   -<colorspace> - option indicating which colorspace color figures should
%                   be saved in: RGB (default), CMYK or gray. CMYK is only
%                   supported in pdf, eps and tiff output.
%   -q<val> - option to vary bitmap image quality (in pdf, eps and jpg
%             files only).  Larger val, in the range 0-100, gives higher
%             quality/lower compression. val > 100 gives lossless
%             compression. Default: '-q95' for jpg, ghostscript prepress
%             default for pdf & eps. Note: lossless compression can
%             sometimes give a smaller file size than the default lossy
%             compression, depending on the type of images.
%   -append - option indicating that if the file (pdfs only) already
%             exists, the figure is to be appended as a new page, instead
%             of being overwritten (default).
%   handle - The handle of the figure or axes to be saved. Default: gcf.
%
%OUT:
%   im - MxNxC uint8 image array of the figure.
%   alpha - MxN single array of alphamatte values in range [0,1], for the
%           case when the background is transparent.
%
%   Some helpful examples and tips can be found at:
%      http://sites.google.com/site/oliverwoodford/software/export_fig
%
%   See also PRINT, SAVEAS.

% Copyright (C) Oliver Woodford 2008-2010

% The idea of using ghostscript is inspired by Peder Axensten's SAVEFIG
% (fex id: 10889) which is itself inspired by EPS2PDF (fex id: 5782).
% The idea for using pdftops came from the MATLAB newsgroup (id: 168171).
% The idea of editing the EPS file to change line styles comes from Jiro
% Doke's FIXPSLINESTYLE (fex id: 17928).
% The idea of changing dash length with line width came from comments on
% fex id: 5743, but the implementation is mine :)
% The idea of anti-aliasing bitmaps came from Anders Brun's MYAA (fex id:
% 20979).
% The idea of appending figures in pdfs came from Matt C in comments on the
% FEX (id: 23629)

% Thanks to Roland Martin for pointing out the colour MATLAB
% bug/feature with colorbar axes and transparent backgrounds.
% Thanks also to Andrew Matthews for describing a bug to do with the figure
% size changing in -nodisplay mode. I couldn't reproduce it, but included a
% fix anyway.
% Thanks to Tammy Threadgill for reporting a bug where an axes is not
% isolated from gui objects.

function [im alpha] = export_fig(varargin)
% Parse the input arguments
[fig options] = parse_args(nargout, varargin{:});
% Isolate the subplot, if it is one
cls = strcmp(get(fig, 'Type'), 'axes');
if cls
    % Given a handle of a single set of axes
    fig = isolate_subplot(fig);
else
    old_mode = get(fig, 'InvertHardcopy');
end
% Hack the font units where necessary (due to a font rendering bug in
% print?). This may not work perfectly in all cases. Also it can change the
% figure layout if reverted, so use a copy.
magnify = options.magnify * options.aa_factor;
if isbitmap(options) && magnify ~= 1
    fontu = findobj(fig, 'FontUnits', 'normalized');
    if ~isempty(fontu)
        % Some normalized font units found
        if ~cls
            fig = copyfig(fig);
            set(fig, 'Visible', 'off');
            fontu = findobj(fig, 'FontUnits', 'normalized');
            cls = true;
        end
        set(fontu, 'FontUnits', 'points');
    end
end
% Set to print exactly what is there
set(fig, 'InvertHardcopy', 'off');
% Set the renderer
switch options.renderer
    case 1
        renderer = '-opengl';
    case 2
        renderer = '-zbuffer';
    case 3
        renderer = '-painters';
    otherwise
        renderer = '-opengl'; % Default for bitmaps
end
% Do the bitmap formats first
if isbitmap(options)
    % Get the background colour
    tcol = get(fig, 'Color');
    if isequal(tcol, 'none') && (options.png || options.alpha)
        % Get out an alpha channel
        % MATLAB "feature": black colorbar axes can change to white and vice versa!
        hCB = findobj(fig, 'Type', 'axes', 'Tag', 'Colorbar');
        if isempty(hCB)
            yCol = [];
            xCol = [];
        else
            yCol = get(hCB, 'YColor');
            xCol = get(hCB, 'XColor');
            if iscell(yCol)
                yCol = cell2mat(yCol);
                xCol = cell2mat(xCol);
            end
            yCol = sum(yCol, 2);
            xCol = sum(xCol, 2);
        end
        % MATLAB "feature": apparently figure size can change when changing
        % colour in -nodisplay mode
        pos = get(fig, 'Position');
        % Set the background colour to black, and set size in case it was
        % changed internally
        set(fig, 'Color', 'k', 'Position', pos);
        % Correct the colorbar axes colours
        set(hCB(yCol==0), 'YColor', [0 0 0]);
        set(hCB(xCol==0), 'XColor', [0 0 0]);
        % Print large version to array
        B = print2array(fig, magnify, renderer);
        % Downscale the image
        B = downsize(single(B), 0, options.aa_factor);
        % Set background to white (and set size)
        set(fig, 'Color', 'w', 'Position', pos);
        % Correct the colorbar axes colours
        set(hCB(yCol==3), 'YColor', [1 1 1]);
        set(hCB(xCol==3), 'XColor', [1 1 1]);
        % Print large version to array
        A = print2array(fig, magnify, renderer);
        % Downscale the image
        A = downsize(single(A), 255, options.aa_factor);
        % Set the background colour (and size) back to normal
        set(fig, 'Color', 'none', 'Position', pos);
        % Compute the alpha map
        alpha = sum(B - A, 3) / (255*3) + 1;
        A = alpha;
        A(A==0) = 1;
        A = B ./ A(:,:,[1 1 1]);
        clear B
        % Convert to greyscale
        if options.colourspace == 2
            A = rgb2grey(A);
        end
        A = uint8(A);
        % Crop the background
        if options.crop
            [alpha v] = crop_background(alpha, 0);
            A = A(v(1):v(2),v(3):v(4),:);
        end
        if options.png
            % Compute the resolution
            res = options.magnify * get(0, 'ScreenPixelsPerInch') / 25.4e-3;
            % Save the png
            imwrite(A, [options.name '.png'], 'Alpha', alpha, 'ResolutionUnit', 'meter', 'XResolution', res, 'YResolution', res);
            % Clear the png bit
            options.png = false;
        end
        % Return only one channel for greyscale
        if isbitmap(options)
            A = check_greyscale(A);
        end
        if options.alpha
            % Store the image
            im = A;
            % Clear the alpha bit
            options.alpha = false;
        end
        % Get the non-alpha image
        if isbitmap(options)
            alph = alpha(:,:,ones(1, size(A, 3)));
            A = uint8(single(A) .* alph + 255 * (1 - alph));
            clear alph
        end
        if options.im
            % Store the new image
            im = A;
        end
    else
        % Print large version to array
        if isequal(tcol, 'none')
            % MATLAB "feature": apparently figure size can change when changing
            % colour in -nodisplay mode
            pos = get(fig, 'Position');
            set(fig, 'Color', 'w', 'Position', pos);
            A = print2array(fig, magnify, renderer);
            set(fig, 'Color', 'none', 'Position', pos);
            tcol = 255;
        else
            A = print2array(fig, magnify, renderer);
            tcol = tcol * 255;
            if ~isequal(tcol, round(tcol))
                tcol = squeeze(A(1,1,:));
            end
        end
        % Crop the background
        if options.crop
            A = crop_background(A, tcol);
        end
        % Downscale the image
        A = downsize(A, tcol, options.aa_factor);
        if options.colourspace == 2
            % Convert to greyscale
            A = rgb2grey(A);
        else
            % Return only one channel for greyscale
            A = check_greyscale(A);
        end
        % Outputs
        if options.im
            im = A;
        end
        if options.alpha
            im = A;
            alpha = zeros(size(A, 1), size(A, 2), 'single');
        end
    end
    % Save the images
    if options.png
        res = options.magnify * get(0, 'ScreenPixelsPerInch') / 25.4e-3;
        imwrite(A, [options.name '.png'], 'ResolutionUnit', 'meter', 'XResolution', res, 'YResolution', res);
    end
    if options.bmp
        imwrite(A, [options.name '.bmp']);
    end
    % Save jpeg with given quality
    if options.jpg
        quality = options.quality;
        if isempty(quality)
            quality = 95;
        end
        if quality > 100
            imwrite(A, [options.name '.jpg'], 'Mode', 'lossless');
        else
            imwrite(A, [options.name '.jpg'], 'Quality', quality);
        end
    end
    % Save tif images in cmyk if wanted (and possible)
    if options.tif
        if options.colourspace == 1 && size(A, 3) == 3
            A = double(255 - A);
            K = min(A, [], 3);
            K_ = 255 ./ max(255 - K, 1);
            C = (A(:,:,1) - K) .* K_;
            M = (A(:,:,2) - K) .* K_;
            Y = (A(:,:,3) - K) .* K_;
            A = uint8(cat(3, C, M, Y, K));
            clear C M Y K K_
        end
        append_mode = {'overwrite', 'append'};
        imwrite(A, [options.name '.tif'], 'Resolution', options.magnify*get(0, 'ScreenPixelsPerInch'), 'WriteMode', append_mode{options.append+1});
    end
end
% Now do the vector formats
if isvector(options)
    % Set the default renderer to painters
    if ~options.renderer
        renderer = '-painters';
    end
    % Generate some filenames
    tmp_nam = [tempname '.eps'];
    if options.pdf
        pdf_nam = [options.name '.pdf'];
    else
        pdf_nam = [tempname '.pdf'];
    end
    % Generate the options for print
    p2eArgs = {renderer};
    if options.colourspace == 1
        p2eArgs = [p2eArgs {'-cmyk'}];
    end
    if ~options.crop
        p2eArgs = [p2eArgs {'-loose'}];
    end
    try
        % Generate an eps
        print2eps(tmp_nam, fig, p2eArgs{:});
        % Generate a pdf
        eps2pdf(tmp_nam, pdf_nam, 1, options.append, options.colourspace==2, options.quality);
    catch
        % Delete the eps
        delete(tmp_nam);
        rethrow(lasterror);
    end
    % Delete the eps
    delete(tmp_nam);
    if options.eps
        try
            % Generate an eps from the pdf
            pdf2eps(pdf_nam, [options.name '.eps']);
        catch
            if ~options.pdf
                % Delete the pdf
                delete(pdf_nam);
            end
            rethrow(lasterror);
        end
        if ~options.pdf
            % Delete the pdf
            delete(pdf_nam);
        end
    end
end
if cls
    % Close the created figure
    close(fig);
else
    % Reset the hardcopy mode
    set(fig, 'InvertHardcopy', old_mode);
end
return
end

function [fig options] = parse_args(nout, varargin)
% Parse the input arguments
% Set the defaults
fig = get(0, 'CurrentFigure');
options = struct('name', 'export_fig_out', ...
                 'crop', true, ...
                 'renderer', 0, ... % 0: default, 1: OpenGL, 2: ZBuffer, 3: Painters
                 'pdf', false, ...
                 'eps', false, ...
                 'png', false, ...
                 'tif', false, ...
                 'jpg', false, ...
                 'bmp', false, ...
                 'colourspace', 0, ... % ): RGB/gray, 1: CMYK, 2: gray
                 'append', false, ...
                 'im', nout == 1, ...
                 'alpha', nout == 2, ...
                 'aa_factor', 4, ...
                 'magnify', 1, ...
                 'quality', []);
native = false; % Set resolution to native of an image

% Go through the other arguments
for a = 1:nargin-1
    if ishandle(varargin{a})
        fig = varargin{a};
    elseif ischar(varargin{a}) && ~isempty(varargin{a})
        if varargin{a}(1) == '-'
            switch lower(varargin{a}(2:end))
                case 'nocrop'
                    options.crop = false;
                case 'opengl'
                    options.renderer = 1;
                case 'zbuffer'
                    options.renderer = 2;
                case 'painters'
                    options.renderer = 3;
                case 'pdf'
                    options.pdf = true;
                case 'eps'
                    options.eps = true;
                case 'png'
                    options.png = true;
                case {'tif', 'tiff'}
                    options.tif = true;
                case {'jpg', 'jpeg'}
                    options.jpg = true;
                case 'bmp'
                    options.bmp = true;
                case 'rgb'
                    options.colourspace = 0;
                case 'cmyk'
                    options.colourspace = 1;
                case {'gray', 'grey'}
                    options.colourspace = 2;
                case {'a1', 'a2', 'a3', 'a4'}
                    options.aa_factor = str2double(varargin{a}(3));
                case 'append'
                    options.append = true;
                case 'native'
                    native = true;
                otherwise
                    val = str2double(regexp(varargin{a}, '(?<=-(m|M|r|R|q|Q))(\d*\.)?\d+(e-?\d+)?', 'match'));
                    if ~isscalar(val)
                        error('option %s not recognised', varargin{a});
                    end
                    switch lower(varargin{a}(2))
                        case 'm'
                            options.magnify = val;
                        case 'r'
                            options.magnify = val ./ get(0, 'ScreenPixelsPerInch');
                        case 'q'
                            options.quality = max(val, 0);
                    end
            end
        else
            name = varargin{a};
            if numel(name) > 3 && name(end-3) == '.' && any(strcmpi(name(end-2:end), {'pdf', 'eps', 'png', 'tif', 'jpg', 'bmp'}))
                options.(lower(name(end-2:end))) = true;
                name = name(1:end-4);
            end
            options.name = name;
        end
    end
end

% Check we have a figure handle
if isempty(fig)
    error('No figure found');
end

% Set the default format
if ~isvector(options) && ~isbitmap(options)
    options.png = true;
end

% If requested, set the resolution to the native vertical resolution of the
% first suitable image found
if native && isbitmap(options)
    % Find a suitable image
    list = findobj(fig, 'Type', 'image', 'Tag', 'export_fig_native');
    if isempty(list)
        list = findobj(fig, 'Type', 'image', 'Visible', 'on');
    end
    for hIm = list(:)'
        % Check height is >= 2
        height = size(get(hIm, 'CData'), 1);
        if height < 2
            continue
        end
        % Account for the image filling only part of the axes, or vice
        % versa
        yl = get(hIm, 'YData');
        if isscalar(yl)
            yl = [yl(1)-0.5 yl(1)+height+0.5];
        else
            if ~diff(yl)
                continue
            end
            yl = yl + [-0.5 0.5] * (diff(yl) / (height - 1));
        end
        hAx = get(hIm, 'Parent');
        yl2 = get(hAx, 'YLim');
        % Find the pixel height of the axes
        oldUnits = get(hAx, 'Units');
        set(hAx, 'Units', 'pixels');
        pos = get(hAx, 'Position');
        set(hAx, 'Units', oldUnits);
        if ~pos(4)
            continue
        end
        % Found a suitable image
        % Account for stretch-to-fill being disabled
        pbar = get(hAx, 'PlotBoxAspectRatio');
        pos = min(pos(4), pbar(2)*pos(3)/pbar(1));
        % Set the magnification to give native resolution
        options.magnify = (height * diff(yl2)) / (pos * diff(yl));
        break
    end
end
return
end

function fh = isolate_subplot(ah, vis)
% Isolate the axes in a figure on their own
% Tag the axes so we can find them in the copy
old_tag = get(ah, 'Tag');
set(ah, 'Tag', 'ObjectToCopy');
% Create a new figure exactly the same as the old one
fh = copyfig(ancestor(ah, 'figure')); %copyobj(ancestor(ah, 'figure'), 0);
if nargin < 2 || ~vis
    set(fh, 'Visible', 'off');
end
% Reset the axes tag
set(ah, 'Tag', old_tag);
% Get all the axes and gui objects
axs = get(fh, 'Children');
% Find the objects to save
ah = findobj(axs, 'Tag', 'ObjectToCopy');
if numel(ah) ~= 1
    close(fh);
    error('Too many axes found');
end
I = true(size(axs));
I(axs==ah) = false;
% Set the axes tag to what it should be
set(ah, 'Tag', old_tag);
% Keep any legends which overlap the subplot
ax_pos = get(ah, 'OuterPosition');
ax_pos(3:4) = ax_pos(3:4) + ax_pos(1:2);
for ah = findobj(axs, 'Tag', 'legend', '-or', 'Tag', 'Colorbar')'
    leg_pos = get(ah, 'OuterPosition');
    leg_pos(3:4) = leg_pos(3:4) + leg_pos(1:2);
    % Overlap test
    if leg_pos(1) < ax_pos(3) && leg_pos(2) < ax_pos(4) &&...
       leg_pos(3) > ax_pos(1) && leg_pos(4) > ax_pos(2)
        I(axs==ah) = false;
    end
end
% Delete all axes except for the input axes and associated items
delete(axs(I));
return
end

function fh = copyfig(fh)
% Is there a legend?
if isempty(findobj(fh, 'Type', 'axes', 'Tag', 'legend'))
    % Safe to copy using copyobj
    fh = copyobj(fh, 0);
else
    % copyobj will change the figure, so save and then load it instead
    tmp_nam = [tempname '.fig'];
    hgsave(fh, tmp_nam);
    fh = hgload(tmp_nam);
    delete(tmp_nam);
end
return
end

function A = downsize(A, padval, factor)
% Downsample an image
if factor == 1
    % Nothing to do
    return
end
try
    % Faster, but requires image processing toolbox
    A = imresize(A, 1/factor, 'bilinear');
catch
    % No image processing toolbox - resize manually
    % Lowpass filter - use Gaussian as is separable, so faster
    switch factor
        case 4
            % sigma: 1.7
            filt = single([0.0148395 0.0498173 0.118323 0.198829 0.236384 0.198829 0.118323 0.0498173 0.0148395]);
        case 3
            % sigma: 1.35
            filt = single([0.025219 0.099418 0.226417 0.297892 0.226417 0.099418 0.025219]);
        case 2
            % sigma: 1.0
            filt = single([0.054489 0.244201 0.40262 0.244201 0.054489]);
    end
    padding = floor(numel(filt) / 2);
    if numel(padval) == 3 && padval(1) == padval(2) && padval(2) == padval(3)
        padval = padval(1);
    end
    if numel(padval) == 1
        B = repmat(single(padval), [size(A, 1) size(A, 2)] + (2 * padding));
    end
    for a = 1:size(A, 3)
        if numel(padval) == 3
            B = repmat(single(padval(a)), [size(A, 1) size(A, 2)] + (2 * padding));
        end
        B(padding+1:end-padding,padding+1:end-padding) = A(:,:,a);
        A(:,:,a) = conv2(filt, filt', B, 'valid');
    end
    clear B
    % Subsample
    A = A(2:factor:end,2:factor:end,:);
end
return
end

function A = rgb2grey(A)
A = cast(reshape(reshape(single(A), [], 3) * single([0.299; 0.587; 0.114]), size(A, 1), size(A, 2)), class(A));
return
end

function A = check_greyscale(A)
% Check if the image is greyscale
if size(A, 3) == 3 && ...
        all(reshape(A(:,:,1) == A(:,:,2), [], 1)) && ...
        all(reshape(A(:,:,2) == A(:,:,3), [], 1))
    A = A(:,:,1); % Save only one channel for 8-bit output
end
return
end

function [A v] = crop_background(A, bcol)
% Map the foreground pixels
[h w c] = size(A);
if isscalar(bcol) && c > 1
    bcol = bcol(ones(1, c));
end
bail = false;
for l = 1:w
    for a = 1:c
        if ~all(A(:,l,a) == bcol(a))
            bail = true;
            break;
        end
    end
    if bail
        break;
    end
end
bail = false;
for r = w:-1:l
    for a = 1:c
        if ~all(A(:,r,a) == bcol(a))
            bail = true;
            break;
        end
    end
    if bail
        break;
    end
end
bail = false;
for t = 1:h
    for a = 1:c
        if ~all(A(t,:,a) == bcol(a))
            bail = true;
            break;
        end
    end
    if bail
        break;
    end
end
bail = false;
for b = h:-1:t
    for a = 1:c
        if ~all(A(b,:,a) == bcol(a))
            bail = true;
            break;
        end
    end
    if bail
        break;
    end
end
% Crop the background
A = A(t:b,l:r,:);
v = [t b l r];
return
end

function b = isvector(options)
b = options.pdf || options.eps;
return
end

function b = isbitmap(options)
b = options.png || options.tif || options.jpg || options.bmp || options.im || options.alpha;
return
end

% Copyright (c) 2009, urs (us) schwarz
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.

% Taken from
% http://www.mathworks.com/matlabcentral/fileexchange/3218-mtit-a-pedestria
% n-major-title-creator

%MTIT		creates a major title in a figure with many axes
%
%		MTIT
%		- creates a major title above all
%		  axes in a figure
%		- preserves the stack order of
%		  the axis handles
%
%SYNTAX
%-------------------------------------------------------------------------------
%		P = MTIT(TXT,[OPT1,...,OPTn])
%		P = MTIT(FH,TXT,[OPT1,...,OPTn])
%
%INPUT
%-------------------------------------------------------------------------------
%    FH	:	a valid figure handle		[def: gcf]
%   TXT	:	title string
%
% OPT	:	argument
% -------------------------------------------
%  xoff	:	+/- displacement along X axis
%  yoff	:	+/- displacement along Y axis
%  zoff	:	+/- displacement along Z axis
%
%		title modifier pair(s)
% -------------------------------------------
%   TPx	:	TVx
%		see: get(text) for possible
%		     parameters/values
%
%OUTPUT
%-------------------------------------------------------------------------------
% par	:	parameter structure
%  .pos :	position of surrounding axis
%   .oh	:	handle of last used axis
%   .ah :	handle of invisible surrounding axis
%   .th :	handle of main title
%
%EXAMPLE
%-------------------------------------------------------------------------------
%	subplot(2,3,[1 3]);		title('PLOT 1');
%	subplot(2,3,4); 		title('PLOT 2');
%	subplot(2,3,5); 		title('PLOT 3');
%	axes('units','inches',...
%	     'color',[0 1 .5],...
%	     'position',[.5 .5 2 2]);	title('PLOT 41');
%	axes('units','inches',...
%	     'color',[0 .5 1],...
%	     'position',[3.5 .5 2 2]);	title('PLOT 42');
%	shg;
%	p=mtit('the BIG title',...
%	     'fontsize',14,'color',[1 0 0],...
%	     'xoff',-.1,'yoff',.025);
% % refine title using its handle <p.th>
%	set(p.th,'edgecolor',.5*[1 1 1]);

% created:
%	us	24-Feb-2003		/ R13
% modified:
%	us	24-Feb-2003		/ CSSM
%	us	06-Apr-2003		/ TMW
%	us	13-Nov-2009 17:38:17

%-------------------------------------------------------------------------------
function	par=mtit(varargin)

		defunit='normalized';
	if	nargout
		par=[];
	end

% check input
	if	nargin < 1
		help(mfilename);
		return;
	end
	if	isempty(get(0,'currentfigure'))
		disp('MTIT> no figure');
		return;
	end

		vl=true(size(varargin));
	if	ischar(varargin{1})
		vl(1)=false;
		figh=gcf;
		txt=varargin{1};
	elseif	any(ishandle(varargin{1}(:)))		&&...
		ischar(varargin{2})
		vl(1:2)=false;
		figh=varargin{1};
		txt=varargin{2};
	else
		error('MTIT> invalid input');
	end
		vin=varargin(vl);
		[off,vout]=get_off(vin{:});

% find surrounding box
		ah=findall(figh,'type','axes');
	if	isempty(ah)
		disp('MTIT> no axis');
		return;
	end
		oah=ah(1);

		ou=get(ah,'units');
		set(ah,'units',defunit);
		ap=get(ah,'position');
	if	iscell(ap)
		ap=cell2mat(get(ah,'position'));
	end
		ap=[	min(ap(:,1)),max(ap(:,1)+ap(:,3)),...
			min(ap(:,2)),max(ap(:,2)+ap(:,4))];
		ap=[	ap(1),ap(3),...
			ap(2)-ap(1),ap(4)-ap(3)];

% create axis...
		xh=axes('position',ap);
% ...and title
		th=title(txt,vout{:});
		tp=get(th,'position');
		set(th,'position',tp+off);
		set(xh,'visible','off','hittest','on');
		set(th,'visible','on');

% reset original units
		ix=find(~strcmpi(ou,defunit));
	if	~isempty(ix)
	for	i=ix(:).'
		set(ah(i),'units',ou{i});
	end
	end

% ...and axis' order
		uistack(xh,'bottom');
		axes(oah);				%#ok

	if	nargout
		par.pos=ap;
		par.oh=oah;
		par.ah=xh;
		par.th=th;
	end
end
%-------------------------------------------------------------------------------
function	[off,vout]=get_off(varargin)

% search for pairs <.off>/<value>

		off=zeros(1,3);
		io=0;
	for	mode={'xoff','yoff','zoff'};
		ix=strcmpi(varargin,mode);
	if	any(ix)
		io=io+1;
		yx=find(ix);
		ix(yx+1)=1;
		off(1,io)=varargin{yx(end)+1};
		varargin=varargin(xor(ix,1));
	end
	end
		vout=varargin;
end
%-------------------------------------------------------------------------------% Copyright (c) 2010, Oliver Woodford
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.
%PDF2EPS  Convert a pdf file to eps format using pdftops
%
% Examples:
%   pdf2eps source dest
%
% This function converts a pdf file to eps format.
%
% This function requires that you have pdftops, from the Xpdf suite of
% functions, installed on your system. This can be downloaded from:
% http://www.foolabs.com/xpdf  
%
%IN:
%   source - filename of the source pdf file to convert. The filename is
%            assumed to already have the extension ".pdf".
%   dest - filename of the destination eps file. The filename is assumed to
%          already have the extension ".eps".

% Copyright (C) Oliver Woodford 2009-2010

% Thanks to Aldebaro Klautau for reporting a bug when saving to
% non-existant directories.

function pdf2eps(source, dest)
% Construct the options string for pdftops
options = ['-q -paper match -eps -level2 "' source '" "' dest '"'];
% Convert to eps using pdftops
[status message] = pdftops(options);
% Check for error
if status
    % Report error
    if isempty(message)
        error('Unable to generate eps. Check destination directory is writable.');
    else
        error(message);
    end
end
% Fix the DSC error created by pdftops
fid = fopen(dest, 'r+');
if fid == -1
    % Cannot open the file
    return
end
fgetl(fid); % Get the first line
str = fgetl(fid); % Get the second line
if strcmp(str(1:min(13, end)), '% Produced by')
    fseek(fid, -numel(str)-1, 'cof');
    fwrite(fid, '%'); % Turn ' ' into '%'
end
fclose(fid);
return
end

% Copyright (c) 2010, Oliver Woodford
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.
function varargout = pdftops(cmd)
%PDFTOPS  Calls a local pdftops executable with the input command
%
% Example:
%   [status result] = pdftops(cmd)
%
% Attempts to locate a pdftops executable, finally asking the user to
% specify the directory pdftops was installed into. The resulting path is
% stored for future reference.
% 
% Once found, the executable is called with the input command string.
%
% This function requires that you have pdftops (from the Xpdf package)
% installed on your system. You can download this from:
% http://www.foolabs.com/xpdf
%
% IN:
%   cmd - Command string to be passed into pdftops.
%
% OUT:
%   status - 0 iff command ran without problem.
%   result - Output from pdftops.

% Copyright: Oliver Woodford, 2009-2010

% Thanks to Jonas Dorn for the fix for the title of the uigetdir window on
% Mac OS.
% Thanks to Christoph Hertel for pointing out a bug in check_xpdf_path
% under linux.

% Call pdftops
[varargout{1:nargout}] = system(sprintf('"%s" %s', xpdf_path, cmd));
return
end

function path = xpdf_path
% Return a valid path
% Start with the currently set path
path = current_xpdf_path;
% Check the path works
if check_xpdf_path(path)
    return
end
% Check whether the binary is on the path
if ispc
    bin = 'pdftops.exe';
else
    bin = 'pdftops';
end
if check_store_xpdf_path(bin)
    path = bin;
    return
end
% Search the obvious places
if ispc
    path = 'C:\Program Files\xpdf\pdftops.exe';
else
    path = '/usr/local/bin/pdftops';
end
if check_store_xpdf_path(path)
    return
end
% Ask the user to enter the path
while 1
    if strncmp(computer,'MAC',3) % Is a Mac
        % Give separate warning as the uigetdir dialogue box doesn't have a
        % title
        uiwait(warndlg('Pdftops not found. Please locate the program, or install xpdf-tools from http://users.phg-online.de/tk/MOSXS/.'))
    end
    base = uigetdir('/', 'Pdftops not found. Please locate the program.');
    if isequal(base, 0)
        % User hit cancel or closed window
        break;
    end
    base = [base filesep];
    bin_dir = {'', ['bin' filesep], ['lib' filesep]};
    for a = 1:numel(bin_dir)
        path = [base bin_dir{a} bin];
        if exist(path, 'file') == 2
            break;
        end
    end
    if check_store_xpdf_path(path)
        return
    end
end
error('pdftops executable not found.');
end

function good = check_store_xpdf_path(path)
% Check the path is valid
good = check_xpdf_path(path);
if ~good
    return
end
% Update the current default path to the path found
if change_value(path, 'current_xpdf_path_str', [mfilename('fullpath') '.m'])
    warning('Path to pdftops executable could not be saved. Enter it manually in pdftops.m.');
    return
end
return
end

function good = check_xpdf_path(path)
% Check the path is valid
[good message] = system(sprintf('"%s" -h', path));
% system returns good = 1 even when the command runs
% Look for something distinct in the help text
good = ~isempty(strfind(message, 'PostScript'));
end

function current_xpdf_path_str = current_xpdf_path
current_xpdf_path_str = 'C:\Program Files\xpdf-3.02pl4-win32\pdftops.exe';
end
% Copyright (c) 2010, Oliver Woodford
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.
%PRINT2ARRAY  Exports a figure to an image array
%
% Examples:
%   A = print2array
%   A = print2array(figure_handle)
%   A = print2array(figure_handle, resolution)
%   A = print2array(figure_handle, resolution, renderer)
%
% This function outputs a bitmap image of the given figure, at the desired
% resolution.
%
% If renderer is '-painters' then ghostcript needs to be installed. This
% can be downloaded from: http://www.ghostscript.com
%
% IN:
%   figure_handle - The handle of the figure to be exported. Default: gcf.
%   resolution - Resolution of the output, as a factor of screen
%                resolution. Default: 1.
%   renderer - string containing the renderer paramater to be passed to
%              print. Default: '-opengl'.
%
% OUT:
%   A - MxNx3 uint8 image of the figure.

% Copyright (C) Oliver Woodford 2008-2010

function A = print2array(fig, res, renderer)
% Generate default input arguments, if needed
if nargin < 2
    res = 1;
    if nargin < 1
        fig = gcf;
    end
end
% Warn if output is large
old_mode = get(fig, 'Units');
set(fig, 'Units', 'pixels');
px = get(fig, 'Position');
set(fig, 'Units', old_mode);
px = prod(px(3:4)*res)/1e6;
if px > 30
    % 30M pixels or larger!
    warning('MATLAB:LargeImage', 'print2array generating a %.1fM pixel image. This could be slow and might also cause memory problems.', px);
end
% Set the resolution parameter
res = ['-r' num2str(ceil(get(0, 'ScreenPixelsPerInch')*res))];
% Generate temporary file name
tmp_nam = [tempname '.tif'];
if nargin > 2 && strcmp(renderer, '-painters')
    % Print to eps file
    tmp_eps = [tempname '.eps'];
    print2eps(tmp_eps, fig, renderer, '-loose');
    try
        % Export to tiff using ghostscript
        ghostscript(['-dEPSCrop -q -dNOPAUSE -dBATCH ' res ' -sDEVICE=tiff24nc -sOutputFile="' tmp_nam '" "' tmp_eps '"']);
    catch
        % Delete the intermediate file
        delete(tmp_eps);
        rethrow(lasterror);
    end
    % Delete the intermediate file
    delete(tmp_eps);
    % Read in the generated bitmap
    A = imread(tmp_nam);
    % Delete the temporary bitmap file
    delete(tmp_nam);
    % Retrieve the background colour
    bcol = get(fig, 'Color');
    % Set border pixels to the correct colour
    if ~isequal(bcol, 'none') && ~isequal(bcol, [1 1 1])
        for l = 1:size(A, 2)
            if ~all(reshape(A(:,l,:) == 255, [], 1))
                break;
            end
        end
        for r = size(A, 2):-1:l
            if ~all(reshape(A(:,r,:) == 255, [], 1))
                break;
            end
        end
        for t = 1:size(A, 1)
            if ~all(reshape(A(t,:,:) == 255, [], 1))
                break;
            end
        end
        for b = size(A, 1):-1:t
            if ~all(reshape(A(b,:,:) == 255, [], 1))
                break;
            end
        end
        bcol = round(bcol * 255);
        for c = 1:size(A, 3)
            A(:,[1:l-1, r+1:end],c) = bcol(c);
            A([1:t-1, b+1:end],:,c) = bcol(c);
        end
    end
else
    if nargin < 3
        renderer = '-opengl';
    end
    % Set paper size
    old_mode = get(fig, 'PaperPositionMode');
    set(fig, 'PaperPositionMode', 'auto');
    try
        % Try hardcopy first - undocumented MATLAB function!
        A = hardcopy(fig, ['-D' renderer(2:end)], res);
    catch
        % Print to tiff file
        print(fig, renderer, res, '-dtiff', tmp_nam);
        % Read in the printed file
        A = imread(tmp_nam);
        % Delete the temporary file
        delete(tmp_nam);
    end
    % Reset paper size
    set(fig, 'PaperPositionMode', old_mode);
end
return
end
% Copyright (c) 2010, Oliver Woodford
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.
%PRINT2EPS  Prints figures to eps with improved line styles
%
% Examples:
%   print2eps filename
%   print2eps(filename, fig_handle)
%   print2eps(filename, fig_handle, options)
%
% This function saves a figure as an eps file, and improves the line style,
% making dashed lines more like those on screen and giving grid lines their
% own dotted style.
%
%IN:
%   filename - string containing the name (optionally including full or
%              relative path) of the file the figure is to be saved as. A
%              ".eps" extension is added if not there already. If a path is
%              not specified, the figure is saved in the current directory.
%   fig_handle - The handle of the figure to be saved. Default: gcf.
%   options - Additional parameter strings to be passed to print.

% Copyright (C) Oliver Woodford 2008-2009

% The idea of editing the EPS file to change line styles comes from Jiro
% Doke's FIXPSLINESTYLE (fex id: 17928)
% The idea of changing dash length with line width came from comments on
% fex id: 5743, but the implementation is mine :)

function print2eps(name, fig, varargin)
options = {'-depsc2'};
if nargin < 2
    fig = gcf;
elseif nargin > 2
    options = [options varargin];
end
% Construct the filename
if numel(name) < 5 || ~strcmpi(name(end-3:end), '.eps')
    name = [name '.eps']; % Add the missing extension
end
% Set paper size
old_mode = get(fig, 'PaperPositionMode');
set(fig, 'PaperPositionMode', 'auto');
% Print to eps file
print(fig, options{:}, name);
% Reset paper size
set(fig, 'PaperPositionMode', old_mode);
% Fix the line styles
try
    fix_lines(name);
catch
    warning('fix_lines failed. This is usually because the figure contains a large number of patch objects. Consider exporting to a bitmap format in this case.');
end
return
end
