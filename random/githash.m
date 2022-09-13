% [hash] = githash(fname, gitdir)

function [hash] = githash(fname, gitdir)

    if ~exist('fname', 'var')
        fname = '.';
    end

    if ~exist('gitdir', 'var')
        gitdir = '';
    else
        gitdir = ['--git-dir=' gitdir];
    end

    if ~findstr(gitdir, '.git')
        gitdir = [gitdir '/.git/']
    end

    [~, hashout] = system(['TERM=xterm git ' gitdir ...
                        ' log -n 1 --no-color --pretty=format:''%H'' ''' fname ''' < /dev/null']);

    % remove bash escape characters
    try
        hash = hashout(9:48)
    catch ME
        warning('githash failed');
        hash = '';
    end
