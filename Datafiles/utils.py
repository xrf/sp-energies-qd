def skip_comment_char(read_func, filename):
    with open(filename) as f:
        s = f.read(2)
        assert s == "# "
        return read_func(f)
