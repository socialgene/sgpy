import tarfile


def untargz(input_path, file_name=False, file_name_ends_with=False):
    """Get a file from a tar.gz
    Parameters
    ----------
    input_path : str
        path to the tar.gz
    file_name : str
        file name
    file_name_ends_with : str
        file ned with
    Returns
    -------
    io.BufferedReader object
    Raises
    -------
    ValueError
        "Too many matches. Function should only pull out a single file"
    """
    # extract domtblout from the tar
    tarfile_object = tarfile.open(input_path)
    temp = tarfile_object.getnames()
    if file_name and file_name_ends_with:
        raise ValueError("Must only set one of 'file_name' or 'file_name_ends_with'")
    if file_name:
        match_result = [x == file_name for x in temp]
        match_result = [i for i, x in enumerate(match_result) if x]
        if len(match_result) > 1:
            raise ValueError(
                "Too many matches. Function should only pull out a single file"
            )
    elif file_name_ends_with:
        # endswith matching because domtblout files in the tar.gz
        # have names that are hashes
        match_result = [x.endswith(file_name_ends_with) for x in temp]
        match_result = [i for i, x in enumerate(match_result) if x]
        if len(match_result) > 1:
            raise ValueError(
                "Too many matches. Function should only pull out a single file"
            )
    else:
        raise ValueError("Must set either 'file_name' or 'file_name_ends_with'")
    file_name_to_extract = temp[match_result[0]]
    file = tarfile_object.extractfile(file_name_to_extract)
    return file
