// Which model is open, and what else can be opened.

import { emit, h, invoke, on, pickFile, pickSave, setModel, solve, state, toast } from '../core.js';

const BIGG = 'Downloads once. BiGG terms: research, education, non-profit.';
const CATALOG = [
  { name: 'iCH360', title: 'E. coli iCH360', detail: '349 reactions. Energy and biosynthesis, on the authors\' map.', terms: 'Built in. Its terms are below.', bundled: true },
  { name: 'e_coli_core', title: 'E. coli core', detail: '95 reactions. The textbook model.', terms: BIGG },
  { name: 'iML1515', title: 'E. coli iML1515', detail: '2,712 reactions. Genome scale.', terms: BIGG },
  { name: 'iJO1366', title: 'E. coli iJO1366', detail: '2,583 reactions. Genome scale.', terms: BIGG },
  { name: 'Human-GEM', title: 'Human-GEM', detail: '12,748 reactions. Human metabolism, reaches the Atlas map.', terms: 'Downloads once. CC BY 4.0.' },
  { name: 'Recon3D', title: 'Recon3D', detail: '10,600 reactions. Human metabolism.', terms: BIGG },
];

let refs = {};

function renderCurrent() {
  const m = state.model;
  if (!m) return;
  const facts = [
    [m.reactions.length, 'Reactions'],
    [m.metabolites.length, 'Metabolites'],
    [m.genes.length, 'Genes'],
    [m.compartments, 'Compartments'],
    [m.balance.balanced, 'Balanced reactions'],
    [m.balance.imbalanced, 'Imbalanced reactions'],
    [m.balance.formula_missing, 'Without formulas'],
  ];
  refs.current.replaceChildren(
    h('h2', {}, m.name || m.id),
    h('div', { class: 'tiles' }, facts.map(([value, label]) => h('div', { class: 'tile' }, h('strong', {}, Number(value).toLocaleString('en')), h('span', {}, label)))),
    h('p', { class: 'hint' }, `Objective: ${m.objective.join(', ') || 'none'}`),
    ...(state.modelPath ? [h('button', { class: 'button button-quiet', type: 'button', onclick: exportSbml }, 'Save as SBML Level 3')] : []),
  );
}

async function exportSbml() {
  const output = await pickSave('Save as SBML', `${state.model.id}.xml`, ['xml']);
  if (!output) return;
  try {
    await invoke('export_sbml', { input: state.modelPath, output });
    toast(`Saved ${output}`);
  } catch (error) {
    toast(String(error));
  }
}

async function openPath(path) {
  refs.status.textContent = 'Opening';
  try {
    const view = await invoke('open_file', { path });
    setModel(view, path);
    await solve();
    refs.status.textContent = '';
    emit('navigate', 'map');
  } catch (error) {
    refs.status.textContent = '';
    toast(String(error));
  }
}

async function openCatalog(item) {
  if (item.bundled) {
    const view = await invoke('open_bundled', { name: item.name });
    setModel(view);
    await solve();
    emit('navigate', 'map');
    return;
  }
  refs.status.textContent = `Getting ${item.title}`;
  try {
    const path = await invoke('download_model', { name: item.name });
    await openPath(path);
  } catch (error) {
    refs.status.textContent = '';
    toast(String(error));
  }
}

export default {
  id: 'models',
  label: 'Models',
  icon: ['M12 3l8 4.5v9L12 21l-8-4.5v-9z', 'M12 12l8-4.5', 'M12 12v9', 'M12 12L4 7.5'],
  mount(el) {
    refs.current = h('div', { class: 'panel current-model' });
    refs.status = h('span', { class: 'hint' });
    el.append(
      h('header', { class: 'view-head' },
        h('h1', {}, 'Models'),
        h('div', { class: 'controls' },
          h('button', { class: 'button button-primary', type: 'button', onclick: async () => { const p = await pickFile('Open an SBML model'); if (p) openPath(p); } }, 'Open SBML file'),
          refs.status)),
      refs.current,
      h('div', { class: 'catalog' }, CATALOG.map((item) => h('button', { class: 'catalog-item', type: 'button', onclick: () => openCatalog(item) },
        h('strong', {}, item.title), h('span', {}, item.detail), item.terms && h('span', { class: 'hint' }, item.terms)))),
    );
    on('model', renderCurrent);
    refs.terms = h('details', { class: 'terms' }, h('summary', {}, 'Terms of the built-in models and data'));
    el.append(refs.terms);
    invoke('bundled_terms').then((terms) => {
      for (const [name, text] of terms) refs.terms.append(h('h3', {}, name), h('pre', {}, text));
    });
  },
};
