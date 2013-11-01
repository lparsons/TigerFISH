function ini_config = parse_ini(varargin)
% Parse ini file and use values to override the default parameters
%
%
%   INI_CONFIG = PARSE_INI(INI_FILENAME, [DEFAULT_INI_FILENAME])
%       Read INI_FILENAME and override default parameters
%       Default paramters are documented in default_parameters.ini
%
%   INPUT
%       INI_FILENAME - ini file readable by the IniConfig class (see class
%       definition for supported syntax).
%
%   OUTPUT
%       INI_CONFIG - Instance of IniConfig class with default parameters
%       except where they are overridden by those specified in INI_FILENAME
%
% ------------------------
% Copyright (c) 2010-2013, Lance R. Parsons <lparsons@princeton.edu>, Nikolai Slavov <nslavov@mit.edu>
% Licensed under the BSD 2-Clause License: http://opensource.org/licenses/BSD-2-Clause
% ------------------------
p = mfilename('fullpath');
[pathstr] = fileparts(p);
default_ini_path = [pathstr filesep '..' filesep 'default_parameters.ini'];

%% Parse Arguments
ip = inputParser;
ip.FunctionName = 'parse_ini';
ip.addOptional('ini_file','',@ischar);
ip.addOptional('default_ini_file',default_ini_path,@(x) exist(x, 'file'));
ip.parse(varargin{:});


%% Default Values
default_ini = IniConfig();
if ~default_ini.ReadFile(ip.Results.default_ini_file)
    err = MException('TigerFISH:Config:ParseError', ...
        sprintf('Unable to read default parameters ini file: "%s"', ip.Results.default_ini_file));
    if ~exist(ip.Results.default_ini_file, 'file')
        errCause = MException('TigerFISH:Config:FileNotFound', ...
            'File not found');
        err = addCause(err, errCause);
    end
    throw(err)
end

%% Override defaults with params specified
ini_config = default_ini;
if ~strcmp(ip.Results.ini_file,'')
    % Custom Values
    custom_ini = IniConfig();
    if custom_ini.ReadFile(ip.Results.ini_file)
        sections = custom_ini.GetSections();
        for s = 1:length(sections)
            [keys, count_keys] = custom_ini.GetKeys(sections{s});
            values = custom_ini.GetValues(sections{s}, keys);
            status = ini_config.SetValues(sections{s}, keys, values);
            if ~status
                err = MException('TigerFISH:Config:ParseError', ...
                    sprintf('Unable to use values from config file: "%s"', ...
                    ip.Results.ini_file));
                throw(err)
            end
        end
    else
        err = MException('TigerFISH:Config:ParseError', ...
            sprintf('Unable to read custom parameters ini file: "%s"', ip.Results.ini_file));
        if ~exist(ip.Results.ini_file, 'file')
            errCause = MException('TigerFISH:Config:FileNotFound', ...
                'File not found');
            err = addCause(err, errCause);
        end
        throw(err)
    end
end
