// The validation results this release was built with.

import { fmt, h, invoke } from '../core.js';

let body;
let loaded = false;

function row(name, reference, result, pass) {
  return h('tr', {},
    h('td', {}, h('span', { class: `dot ${pass ? 'dot-pass' : 'dot-fail'}`, role: 'img', 'aria-label': pass ? 'Passes' : 'Fails' }), name),
    h('td', {}, reference),
    h('td', {}, result));
}

const count = (v) => Number(v).toLocaleString('en');

function rows(v) {
  const core = v.v01_core_cobrapy.checks;
  const keio = v.v02_iml1515_keio;
  const human = v.v03_human_gem;
  const teusink = v.v04_teusink_mca;
  const suite = v.v05_sbml_suite;
  const sampling = v.v06_sampling;
  const overlay = v.v07_overlay_escher;
  const numerics = v.v08_numerics;
  const recon = v.v11_recon3d.checks;
  const enzyme = v.v12_ecmodel;
  const worst = Math.max(...Object.values(teusink.theorems).map((t) => t.worst));
  const failingCases = suite.stochastic_cases_failing_the_literal_criterion ?? [];
  return [
    row('Growth of e_coli_core', 'BiGG 0.8739 1/h, COBRApy', core.fba_objective.rashevsky_highs.toFixed(6), core.fba_objective.pass),
    row('Flux variability', 'COBRApy within 1e-6', `Largest difference ${fmt.sci(core['fva_fraction_1.0'].max_abs_difference)}`, core['fva_fraction_1.0'].pass),
    row('Double gene deletions', 'COBRApy within 1e-5', `${count(core.double_gene_deletion.pairs)} pairs, largest difference ${fmt.sci(core.double_gene_deletion.max_abs_difference)}`, core.double_gene_deletion.pass),
    row('iML1515 essentiality', 'Keio screens, Monk et al. 2017, 93.4 %', `${keio.accuracy_percent.toFixed(2)} % of ${count(keio.phenotypes)} phenotypes`, keio.pass),
    row('Human-GEM tasks', 'All 57 essential tasks, Robinson et al. 2020', `${human.essential_passed} of ${human.essential_tasks}, same verdicts as RAVEN`, human.pass),
    row('Recon3D', 'COBRApy, 2,248 gene deletions within 1e-6', `Objective ${recon.objective.rashevsky.toFixed(4)}, largest difference ${fmt.sci(recon.gene_deletions.max_relative_difference)}`, v.v11_recon3d.pass),
    row('Enzyme constraints', 'Published eciML1515, certified optimum in 7 conditions', `Largest difference ${fmt.sci(Math.max(...enzyme.conditions.map((c) => c.module_relative_difference)))}`, enzyme.pass),
    row('Yeast glycolysis', 'Teusink et al. 2000 Table 4, libRoadRunner', `Largest relative difference ${fmt.sci(teusink.steady_state_max_rel_difference_species)}`, teusink.pass),
    row('Control coefficients', 'COPASI within 1e-4', `Largest difference ${fmt.sci(teusink.flux_control_max_abs_difference_vs_copasi)}`, teusink.pass),
    row('Control theorems', 'Hold to 1e-8 on five models', `Worst residual ${fmt.sci(worst)}`, worst < 1e-8),
    row('SBML semantic suite, Rodas4', 'At least 95 % of applicable cases', `${count(suite.rodas4.passed)} of ${count(suite.rodas4.applicable)}`, suite.rodas4.percent >= 95),
    row('SBML semantic suite, BDF', 'At least 95 % of applicable cases', `${count(suite.bdf.passed)} of ${count(suite.bdf.applicable)}`, suite.bdf.percent >= 95),
    row('Stochastic suite', 'DSMTS mean and variance criteria', failingCases.length ? `${suite.stochastic_cases - failingCases.length} of ${suite.stochastic_cases}; case ${failingCases.join(', ')} checked against its exact law` : `All ${suite.stochastic_cases}`, suite.pass),
    row('Flux sampling', 'COBRApy OptGP means within 0.15 SD', `Largest difference ${sampling.max_mean_difference_in_flux_sd.chrr_vs_cobrapy.value.toFixed(3)} SD, R-hat ${sampling.chrr_diagnostics.max_split_rhat.toFixed(3)}`, sampling.pass),
    row('Thermodynamics', 'Flux follows the sign of the Gibbs energy', `${numerics.tfa.reactions_with_energies} reactions, none against`, numerics.tfa.pass),
    row('Atlas overlay', 'Independent recomputation, 300 exact events', `${count(overlay.coverage.events_exact)} exact events`, overlay.pass),
    row('Escher output', 'Escher 1-0-0 schema', overlay.escher_schema_valid ? 'Valid' : 'Invalid', overlay.escher_schema_valid),
  ];
}

export default {
  id: 'checks',
  label: 'Checks',
  icon: ['M20 6L9 17l-5-5'],
  mount(el) {
    body = h('tbody');
    el.append(
      h('header', { class: 'view-head' }, h('h1', {}, 'Checks')),
      h('div', { class: 'panel table-panel checks-panel' },
        h('table', { class: 'table' },
          h('thead', {}, h('tr', {}, h('th', { scope: 'col' }, 'What is checked'), h('th', { scope: 'col' }, 'Against'), h('th', { scope: 'col' }, 'Result'))),
          body)),
    );
  },
  async show() {
    if (loaded) return;
    loaded = true;
    body.replaceChildren(...rows(await invoke('validation')));
  },
};
