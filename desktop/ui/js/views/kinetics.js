// Kinetic models: time courses, steady states and control coefficients, stochastic ensembles.

import { fmt, h, invoke, pickFile, state, toast } from '../core.js';
import { heatmap, lineChart } from '../chart.js';

const COLOURS = ['#f4d64e', '#6d97ff', '#7fd6c2', '#ef8a5a', '#c79bff', '#e7e4dc', '#9ad05f', '#ff9ac2', '#58c7f3', '#d8b26a', '#b0b8c4', '#ff6f6f'];
let refs = {};
let opened = false;
let course = null;
let control = null;
let chosen = new Set();
let tab = 'course';
const overrides = new Map();
let hoverHeat = null;

async function open(path = null) {
  try {
    state.kinetic = await invoke('kinetic_open', { path });
    overrides.clear();
    const dynamic = state.kinetic.species.filter((s) => s.dynamic).map((s) => s.id);
    chosen = new Set(dynamic.slice(0, 7));
    refs.title.textContent = state.kinetic.name || state.kinetic.id;
    renderSpecies();
    renderParameters();
    opened = true;
    course = null;
    control = null;
    await simulate();
  } catch (error) {
    toast(String(error));
  }
}

function colourOf(id) {
  const all = state.kinetic.species.map((s) => s.id);
  return COLOURS[all.indexOf(id) % COLOURS.length];
}

function renderSpecies() {
  refs.species.replaceChildren(...state.kinetic.species.map((s) => {
    const on = chosen.has(s.id);
    return h('button', {
      class: `species ${on ? 'is-on' : ''} ${s.dynamic ? '' : 'is-fixed'}`,
      type: 'button',
      'aria-pressed': String(on),
      onclick: () => {
        if (chosen.has(s.id)) chosen.delete(s.id);
        else chosen.add(s.id);
        renderSpecies();
        drawCourse();
      },
    }, h('span', { class: 'species-dot', style: { background: on ? colourOf(s.id) : 'transparent' } }), s.id);
  }));
}

function renderParameters() {
  const q = refs.paramFilter.value.trim().toLowerCase();
  const list = state.kinetic.parameters.filter((p) => !q || p.name.toLowerCase().includes(q)).slice(0, 80);
  refs.params.replaceChildren(...list.map((p) => h('label', { class: 'param' },
    h('span', { class: 'ident' }, p.name),
    h('input', {
      class: 'number',
      type: 'number',
      step: 'any',
      value: overrides.get(p.name) ?? p.value,
      onchange: (e) => {
        overrides.set(p.name, Number(e.target.value));
        rerun();
      },
    }))));
}

async function rerun() {
  if (tab === 'course') await simulate();
  else if (tab === 'control') await runControl();
  else await runStochastic();
}

async function simulate() {
  try {
    refs.status.textContent = 'Integrating';
    course = await invoke('kinetic_simulate', {
      end: Number(refs.end.value),
      steps: 400,
      method: refs.method.value,
      overrides: [...overrides.entries()],
    });
    refs.status.textContent = `${course.stats.steps} steps, ${course.milliseconds.toFixed(0)} ms`;
    drawCourse();
  } catch (error) {
    refs.status.textContent = '';
    toast(String(error));
  }
}

function drawCourse() {
  if (!course) return;
  const series = [...chosen].filter((id) => course.series[id]).map((id) => ({ x: course.time, y: course.series[id], colour: colourOf(id), name: id }));
  lineChart(refs.coursePlot, { series, xlabel: 'Time', ylabel: 'Concentration' });
}

async function runControl() {
  try {
    refs.status.textContent = 'Solving the steady state';
    control = await invoke('kinetic_control', { overrides: [...overrides.entries()] });
    refs.status.textContent = `Residual ${fmt.sci(control.residual)}`;
    renderControl();
  } catch (error) {
    refs.status.textContent = '';
    refs.controlFacts.replaceChildren(h('p', { class: 'hint' }, String(error)));
  }
}

function renderControl() {
  if (!control) return;
  const t = control.theorems;
  const rowsAre = refs.matrix.value;
  const values = rowsAre === 'flux' ? control.flux_control : control.concentration_control;
  const rows = rowsAre === 'flux' ? control.reactions : control.states;
  hoverHeat = heatmap(refs.heat, { values, rows, cols: control.reactions, limit: 1, left: 78, top: 72 });
  refs.controlFacts.replaceChildren(
    h('div', { class: 'tile' }, h('strong', {}, fmt.sci(t.flux_summation)), h('span', {}, 'Flux summation')),
    h('div', { class: 'tile' }, h('strong', {}, fmt.sci(t.concentration_summation)), h('span', {}, 'Concentration summation')),
    h('div', { class: 'tile' }, h('strong', {}, fmt.sci(t.flux_connectivity)), h('span', {}, 'Flux connectivity')),
    h('div', { class: 'tile' }, h('strong', {}, fmt.sci(t.concentration_connectivity)), h('span', {}, 'Concentration connectivity')),
    h('div', { class: 'tile' }, h('strong', {}, String(control.conserved_moieties)), h('span', {}, 'Conserved moieties')),
  );
}

async function runStochastic() {
  try {
    refs.status.textContent = 'Simulating the ensemble';
    const run = await invoke('kinetic_stochastic', {
      end: Number(refs.end.value),
      steps: 100,
      runs: Number(refs.runs.value),
      seed: 1,
      epsilon: refs.leap.checked ? 0.03 : null,
      overrides: [...overrides.entries()],
    });
    refs.status.textContent = `${run.runs.toLocaleString('en')} realizations`;
    const ids = [...chosen].filter((id) => run.mean[id]);
    lineChart(refs.stochPlot, {
      series: ids.map((id) => ({ x: run.time, y: run.mean[id], colour: colourOf(id), name: id })),
      bands: ids.map((id) => ({
        x: run.time,
        low: run.mean[id].map((v, k) => v - run.sd[id][k]),
        high: run.mean[id].map((v, k) => v + run.sd[id][k]),
        colour: colourOf(id),
      })),
      xlabel: 'Time',
      ylabel: 'Amount, mean and one standard deviation',
    });
  } catch (error) {
    refs.status.textContent = '';
    toast(String(error));
  }
}

function setTab(next) {
  tab = next;
  for (const [id, node] of Object.entries(refs.panels)) node.hidden = id !== next;
  refs.tabs.querySelectorAll('[role="tab"]').forEach((b) => b.setAttribute('aria-selected', String(b.dataset.tab === next)));
  requestAnimationFrame(rerun);
}

export default {
  id: 'kinetics',
  label: 'Kinetics',
  icon: ['M4 20V4', 'M4 20h16', 'M4 16c4 0 4-10 8-10s4 12 8 12'],
  mount(el) {
    refs.title = h('h2', { class: 'kinetic-title' }, '');
    refs.end = h('input', { class: 'number', type: 'number', min: 0.001, step: 'any', value: 20, 'aria-label': 'End time', onchange: () => rerun() });
    refs.method = h('select', { class: 'select', 'aria-label': 'Integrator', onchange: () => rerun() }, h('option', { value: 'rodas4' }, 'Rodas4'), h('option', { value: 'bdf' }, 'BDF'));
    refs.status = h('span', { class: 'hint' });
    refs.species = h('div', { class: 'species-list' });
    refs.paramFilter = h('input', { class: 'find', type: 'search', placeholder: 'Filter parameters', 'aria-label': 'Filter parameters', oninput: renderParameters });
    refs.params = h('div', { class: 'param-list' });
    refs.coursePlot = h('canvas', { class: 'plot', 'aria-label': 'Time course' });
    refs.stochPlot = h('canvas', { class: 'plot', 'aria-label': 'Stochastic ensemble' });
    refs.heat = h('canvas', { class: 'plot heat', 'aria-label': 'Control coefficients' });
    refs.heatTip = h('div', { class: 'tooltip', hidden: true });
    refs.heatTip.hidden = true;
    refs.matrix = h('select', { class: 'select', 'aria-label': 'Which coefficients', onchange: renderControl }, h('option', { value: 'flux' }, 'Flux control'), h('option', { value: 'concentration' }, 'Concentration control'));
    refs.controlFacts = h('div', { class: 'tiles' });
    refs.runs = h('select', { class: 'select', 'aria-label': 'Realizations', onchange: () => rerun() }, ...[100, 1000, 10000].map((n) => h('option', { value: n, selected: n === 100 }, n.toLocaleString('en'))));
    refs.leap = h('input', { type: 'checkbox', onchange: () => rerun(), 'aria-label': 'Tau-leaping' });
    refs.panels = {
      course: h('div', { class: 'tab-panel plot-panel panel' }, refs.coursePlot),
      control: h('div', { class: 'tab-panel control-panel', hidden: true },
        h('div', { class: 'panel-head' }, refs.matrix, refs.controlFacts),
        h('div', { class: 'panel plot-panel heat-wrap' }, refs.heat, refs.heatTip)),
      stochastic: h('div', { class: 'tab-panel', hidden: true },
        h('div', { class: 'panel-head' }, h('label', { class: 'field field-inline' }, 'Realizations', refs.runs), h('label', { class: 'field field-inline' }, refs.leap, 'Tau-leaping')),
        h('div', { class: 'panel plot-panel' }, refs.stochPlot)),
    };
    refs.panels.control.hidden = true;
    refs.panels.stochastic.hidden = true;
    refs.tabs = h('div', { class: 'tabs', role: 'tablist' },
      [['course', 'Time course'], ['control', 'Steady state and control'], ['stochastic', 'Stochastic']].map(([id, label]) =>
        h('button', { role: 'tab', type: 'button', class: 'tab', dataset: { tab: id }, 'aria-selected': String(id === 'course'), onclick: () => setTab(id) }, label)));
    el.append(
      h('header', { class: 'view-head' },
        h('h1', {}, 'Kinetics'),
        refs.tabs,
        h('div', { class: 'controls' },
          h('label', { class: 'field' }, 'End time', refs.end),
          h('label', { class: 'field' }, 'Integrator', refs.method),
          h('button', { class: 'button button-quiet', type: 'button', onclick: async () => { const p = await pickFile('Open a kinetic model'); if (p) open(p); } }, 'Open SBML'),
          refs.status)),
      h('div', { class: 'kinetic-layout' },
        h('aside', { class: 'kinetic-side' }, refs.title, h('h3', {}, 'Species'), refs.species, h('h3', {}, 'Parameters'), refs.paramFilter, refs.params),
        h('div', { class: 'kinetic-main' }, ...Object.values(refs.panels))),
    );
    new ResizeObserver(() => { if (tab === 'course') drawCourse(); else if (tab === 'control') renderControl(); }).observe(refs.coursePlot);
    refs.heat.addEventListener('pointermove', (e) => {
      if (!hoverHeat) return;
      const box = refs.heat.getBoundingClientRect();
      const hit = hoverHeat(e.clientX - box.left, e.clientY - box.top);
      refs.heatTip.hidden = !hit;
      if (hit) {
        refs.heatTip.style.left = `${e.clientX - box.left + 12}px`;
        refs.heatTip.style.top = `${e.clientY - box.top + 12}px`;
        refs.heatTip.replaceChildren(h('strong', {}, `${hit.col} on ${hit.row}`), h('span', {}, hit.value === null ? 'Undefined at zero flux' : hit.value.toFixed(4)));
      }
    });
    refs.heat.addEventListener('pointerleave', () => { refs.heatTip.hidden = true; });
  },
  show() {
    if (!opened) open();
    else requestAnimationFrame(rerun);
  },
};
