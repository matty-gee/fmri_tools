function on_path = check_path(check)
    on_path = contains([pathsep, path, pathsep], [pathsep, check, pathsep], 'IgnoreCase', ispc);