def clean_up(ggmv):
    # remove all objects
    for name, entity in ggmv._artists_dict.items():
        ggmv._artists_dict[name] = []