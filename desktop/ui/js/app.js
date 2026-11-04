// Application shell: the rail of workspaces and the opening of the first model.

import { emit, h, icon, invoke, on, setModel, solve, state, toast } from './core.js';
import mapView from './views/map.js';
import screenView from './views/screen.js';
import rangesView from './views/ranges.js';
import samplesView from './views/samples.js';
import kineticsView from './views/kinetics.js';
import atlasView from './views/atlas.js';
import modelsView from './views/models.js';
import checksView from './views/checks.js';

const views = [mapView, screenView, rangesView, samplesView, kineticsView, atlasView, modelsView, checksView];
let active = null;
const buttons = new Map();

function show(id) {
  const view = views.find((v) => v.id === id) ?? views[0];
  if (active === view) return;
  if (active) {
    active.node.hidden = true;
    active.hide?.();
    buttons.get(active.id).removeAttribute('aria-current');
  }
  active = view;
  view.node.hidden = false;
  buttons.get(view.id).setAttribute('aria-current', 'page');
  view.show?.();
  emit('view', view.id);
}

function build() {
  const root = document.getElementById('app');
  const rail = h('nav', { class: 'rail', 'aria-label': 'Workspaces' },
    h('div', { class: 'brand', title: 'Rashevsky' }, h('img', { src: 'icon.svg', alt: 'Rashevsky', width: 40, height: 40 })));
  const stage = h('main', { class: 'stage' });
  views.forEach((view, index) => {
    const button = h('button', { class: 'rail-item', type: 'button', title: `${view.label} (${index + 1})`, onclick: () => show(view.id) },
      icon(view.icon), h('span', {}, view.label));
    buttons.set(view.id, button);
    rail.append(button);
    view.node = h('section', { class: `view view-${view.id}`, hidden: true, 'aria-label': view.label });
    view.node.hidden = true;
    stage.append(view.node);
    view.mount(view.node);
  });
  root.append(rail, stage, h('div', { id: 'toast', class: 'toast', role: 'status', hidden: true }));
  document.getElementById('toast').hidden = true;
  window.addEventListener('keydown', (e) => {
    if (e.target.matches('input, select, textarea') || e.ctrlKey || e.metaKey || e.altKey) return;
    const n = Number(e.key);
    if (n >= 1 && n <= views.length) show(views[n - 1].id);
  });
  on('navigate', show);
}

async function start() {
  build();
  state.startup = await invoke('startup');
  if (state.startup.qa) {
    window.addEventListener('error', (e) => invoke('qa_log', { line: `error ${e.message} at ${e.filename}:${e.lineno}` }));
    window.addEventListener('unhandledrejection', (e) => invoke('qa_log', { line: `rejection ${e.reason}` }));
  }
  show(state.startup.view ?? 'map');
  try {
    if (state.startup.file) {
      const view = await invoke('open_file', { path: state.startup.file });
      setModel(view);
    } else {
      const view = await invoke('open_bundled', { name: 'iCH360' });
      setModel(view);
    }
    await solve();
  } catch (error) {
    toast(`Could not open the model: ${error}`);
    const view = await invoke('open_bundled', { name: 'iCH360' });
    setModel(view);
    await solve();
  }
  // Canvas text uses the font as soon as it has loaded.
  document.fonts?.ready.then(() => emit('fonts'));
  window.__rashevsky = { state, show, solve, emit };
}

start();
