import BioSimSpace as BSS

# Override the equilibration protocol with some custom defaults.


def EquilibrationNVT(
    runtime=5 * BSS.Units.Time.nanosecond,
    report_interval=2500,
    restart_interval=250000,
    restraint="backbone",
):
    return BSS.Protocol.Equilibration(
        runtime=runtime,
        report_interval=report_interval,
        restart_interval=restart_interval,
        restraint=restraint,
    )


def EquilibrationNPT(
    runtime=5 * BSS.Units.Time.nanosecond,
    pressure=BSS.Units.Pressure.atm,
    report_interval=2500,
    restart_interval=250000,
    restraint="backbone",
):
    return BSS.Protocol.Equilibration(
        runtime=runtime,
        pressure=pressure,
        report_interval=report_interval,
        restart_interval=restart_interval,
        restraint=restraint,
    )
