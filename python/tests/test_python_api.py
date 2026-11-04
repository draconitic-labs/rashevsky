"""End-to-end checks of the Python interface on the bundled model."""
import pathlib

import pytest

import rashevsky as rv

CORE = str(pathlib.Path(__file__).resolve().parents[2] / "models" / "e_coli_core" / "e_coli_core.xml.gz")
BIOMASS = "BIOMASS_Ecoli_core_w_GAM"


@pytest.fixture()
def model():
    return rv.read_sbml(CORE)


def test_fba_growth(model):
    assert model.fba()["objective"] == pytest.approx(0.8739215, abs=1e-6)


def test_knockout_of_essential_reaction(model):
    # Without enolase glucose cannot reach pyruvate, which is lethal on this medium.
    assert model.fba(knockout_reactions=["ENO"])["objective"] == pytest.approx(0.0, abs=1e-7)


def test_copy_is_independent(model):
    other = model.copy()
    other.set_bounds("EX_glc__D_e", 0.0, 1000.0)
    assert model.bounds("EX_glc__D_e") == (-10.0, 1000.0)
    # Without glucose the maintenance requirement cannot be met, so the model has no steady state.
    assert other.fba()["status"] == "Infeasible"
    assert model.fba()["status"] == "Optimal"


def test_fva_holds_growth(model):
    ranges = model.fva(fraction=1.0, reactions=[BIOMASS])
    assert ranges[BIOMASS][0] == pytest.approx(0.8739215, abs=1e-6)
    assert ranges[BIOMASS][1] == pytest.approx(0.8739215, abs=1e-6)


def test_pfba_keeps_growth_and_reduces_flux(model):
    plain = model.fba()
    parsimonious = model.pfba()
    assert parsimonious["fluxes"][BIOMASS] == pytest.approx(plain["objective"], abs=1e-6)
    assert sum(abs(v) for v in parsimonious["fluxes"].values()) <= sum(abs(v) for v in plain["fluxes"].values()) + 1e-6


def test_gene_deletion_counts(model):
    assert len(model.single_gene_deletion()) == 137


def test_envelope_and_medium(model):
    points = model.production_envelope("EX_o2_e", "EX_ac_e", points=5)
    assert len(points) == 5
    medium = model.medium()
    assert medium["EX_glc__D_e"] == 10.0
    assert model.minimal_medium(0.1)["total"] > 0


def test_unknown_identifier_raises(model):
    with pytest.raises(ValueError):
        model.set_bounds("no_such_reaction", 0.0, 1.0)


def test_sampling_returns_diagnostics(model):
    result = model.sample(chains=2, samples=100, thinning=20, burn_in=200, fraction=0.9)
    assert result["dimension"] > 0
    assert result["diagnostics"]["max_split_rhat"] < 2.0


SIMPLE = """<?xml version="1.0"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">
<model id="decay"><listOfCompartments><compartment id="c" size="1" constant="true"/></listOfCompartments>
<listOfSpecies><species id="A" compartment="c" initialAmount="1" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/></listOfSpecies>
<listOfParameters><parameter id="k" value="0.5" constant="true"/></listOfParameters>
<listOfReactions><reaction id="r" reversible="false"><listOfReactants><speciesReference species="A" stoichiometry="1" constant="true"/></listOfReactants>
<kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML"><apply><times/><ci>k</ci><ci>A</ci><ci>c</ci></apply></math></kineticLaw></reaction></listOfReactions>
</model></sbml>"""


@pytest.fixture()
def decay(tmp_path):
    path = tmp_path / "decay.xml"
    path.write_text(SIMPLE)
    return rv.KineticModel.read_sbml(str(path))


@pytest.mark.parametrize("method", ["rodas4", "bdf"])
def test_simulation_matches_exact_decay(decay, method):
    import math

    run = decay.simulate(end=2.0, steps=4, method=method, symbols=["A"])
    for t, value in zip(run["time"], run["series"]["A"]):
        assert value == pytest.approx(math.exp(-0.5 * t), rel=1e-6)


def test_fit_recovers_rate(decay):
    import math

    data = [{"time": t, "symbol": "A", "value": math.exp(-0.5 * t), "sigma": 0.01} for t in (0.5, 1.0, 2.0, 4.0)]
    fitted = decay.fit([{"name": "k", "start": 2.0}], data)
    assert fitted["estimates"][0] == pytest.approx(0.5, rel=1e-4)


def test_stochastic_mean(decay, tmp_path):
    path = tmp_path / "many.xml"
    path.write_text(SIMPLE.replace('initialAmount="1"', 'initialAmount="100"'))
    model = rv.KineticModel.read_sbml(str(path))
    result = model.ssa(end=1.0, steps=1, runs=4000, seed=5)
    import math

    assert result["mean"]["A#a"][1] == pytest.approx(100 * math.exp(-0.5), abs=0.5)


BIRTH_DEATH = """<?xml version="1.0"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">
<model id="birth_death"><listOfCompartments><compartment id="c" size="1" constant="true"/></listOfCompartments>
<listOfSpecies><species id="S" compartment="c" initialConcentration="1" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/></listOfSpecies>
<listOfParameters><parameter id="k0" value="1.5" constant="true"/><parameter id="k1" value="2" constant="true"/></listOfParameters>
<listOfReactions>
<reaction id="birth" reversible="false"><listOfProducts><speciesReference species="S" stoichiometry="1" constant="true"/></listOfProducts>
<kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML"><apply><times/><ci>k0</ci><ci>c</ci></apply></math></kineticLaw></reaction>
<reaction id="death" reversible="false"><listOfReactants><speciesReference species="S" stoichiometry="1" constant="true"/></listOfReactants>
<kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML"><apply><times/><ci>k1</ci><ci>S</ci><ci>c</ci></apply></math></kineticLaw></reaction>
</listOfReactions></model></sbml>"""


@pytest.fixture()
def birth_death(tmp_path):
    path = tmp_path / "birth_death.xml"
    path.write_text(BIRTH_DEATH)
    return rv.KineticModel.read_sbml(str(path))


def test_local_sensitivity_equals_the_derivative(birth_death):
    local = birth_death.local_sensitivity(["k0", "k1"])
    assert local["state"][0][0] == pytest.approx(0.5, rel=1e-8)
    assert local["state"][0][1] == pytest.approx(-0.375, rel=1e-8)


def test_sobol_ranks_the_wider_range_first(birth_death):
    ranges = [{"name": "k0", "low": 1.4, "high": 1.6}, {"name": "k1", "low": 1.0, "high": 4.0, "log": True}]
    result = birth_death.sobol(ranges, "S", samples=400, seed=3)
    assert result["failures"] == 0
    assert result["total"][1] > 0.8 > result["total"][0]


def test_sobol_rejects_an_unknown_output(birth_death):
    with pytest.raises(ValueError):
        birth_death.sobol([{"name": "k0", "low": 1.0, "high": 2.0}], "nope", samples=8)
