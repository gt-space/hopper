classdef CodegenTextTest < matlab.unittest.TestCase
    % Code that MATLAB Function blocks compile keeps string literals ASCII.
    % Code generation re-encodes strings as UTF-8; on Linux a character
    % such as '×' then changes length and the block fails to parse (CI
    % failed this way on qmult.m's assert message, which Windows accepted).
    % Comments may use any character.

    properties (Constant)
        % folders whose functions run inside MATLAB Function blocks
        Folders = {'subsystems', fullfile('src', 'navigation'), fullfile('src', 'sensors'), ...
            fullfile('src', 'fdir'), fullfile('src', 'faults'), fullfile('src', 'gnc'), ...
            fullfile('src', 'dynamics'), fullfile('src', 'slosh'), fullfile('src', 'propulsion')}
    end

    methods (Test)
        function stringLiteralsAreAscii(tc)
            bad = {};
            for d = CodegenTextTest.Folders
                files = dir(fullfile(hopperRoot(d{1}), '**', '*.m'));
                for f = files'
                    txt = fileread(fullfile(f.folder, f.name));
                    lines = splitlines(txt);
                    for k = 1:numel(lines)
                        if CodegenTextTest.hasNonAsciiString(lines{k})
                            bad{end+1} = sprintf('%s:%d', erase(fullfile(f.folder, f.name), [hopperRoot() filesep]), k); %#ok<AGROW>
                        end
                    end
                end
            end
            tc.verifyEmpty(bad, sprintf('Non-ASCII characters in string literals:\n  %s', strjoin(bad, '\n  ')));
        end
    end

    methods (Static)
        function tf = hasNonAsciiString(line)
            % a quoted string, before any comment, holding a character > 127
            tf = false;
            inStr = ''; prev = ' ';
            for c = line
                if isempty(inStr)
                    if c == '%', return, end
                    % ' starts a string unless it transposes (after a name, ) ] } or ')
                    if c == '"' || (c == '''' && ~(isstrprop(prev, 'alphanum') || any(prev == ')]}._''')))
                        inStr = c;
                    end
                elseif c == inStr
                    inStr = '';
                elseif double(c) > 127
                    tf = true; return
                end
                prev = c;
            end
        end
    end
end
