import BioSimSpace as BSS


# Wrap the instantiation of BSS.Process.Amber objects to configure them
# such that coordinates are always wrapped to the minimum image.
def Amber(
    system, protocol, exe=None, name="amber", work_dir=None, seed=None, property_map={}
):
    # Create process using the passed parameters.
    process = BSS.Process.Amber(
        system,
        protocol,
        exe=exe,
        name=name,
        work_dir=work_dir,
        seed=seed,
        property_map=property_map,
    )

    # Get the config.
    config = process.getConfig()

    # Add coordinate wrapping as the second last line.
    config[-1] = "  iwrap=1,"
    config.append(" /")

    # Set the new config.
    process.setConfig(config)

    # Return the process.
    return process
