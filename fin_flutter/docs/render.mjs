// Renders docs/figures/src/*.mmd (Mermaid) to docs/figures/*.svg with Pele.
//     node docs/render.mjs [name ...]
// Pele comes from the guide's node_modules (cd flightModel/docs/guide && npm install),
// or from the directory in $PELE_DIR.
import { mkdirSync, readFileSync, readdirSync, writeFileSync } from 'node:fs';
import { basename, dirname, join, resolve } from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';

const here = dirname(fileURLToPath(import.meta.url));
const peleDir = process.env.PELE_DIR ||
  resolve(here, '../../../CFD-ITBA-ROCKETRY/flightModel/docs/guide/node_modules/pele');
const candidates = [peleDir, resolve(here, '../../flightModel/docs/guide/node_modules/pele')];
let render;
for (const dir of candidates) {
  try {
    const pkg = JSON.parse(readFileSync(join(dir, 'package.json'), 'utf8'));
    ({ render } = await import(pathToFileURL(join(dir, pkg.module || pkg.main)).href));
    break;
  } catch { /* next */ }
}
if (!render) {
  console.error('no pele: run npm install in flightModel/docs/guide, or set PELE_DIR');
  process.exit(1);
}

const theme = '--pele-bg:#ffffff;--pele-fg:#1a1a1a;--pele-muted:#5b6573;--pele-line:#1f3a5f;' +
  '--pele-surface:#eaf2fb;--pele-surface-alt:#fff4e0;--pele-border:#2f6db5;--pele-accent:#c77700;';
const src = join(here, 'figures/src');
const out = join(here, 'figures');
mkdirSync(out, { recursive: true });
const names = process.argv.length > 2
  ? process.argv.slice(2)
  : readdirSync(src).filter((f) => f.endsWith('.mmd')).map((f) => basename(f, '.mmd'));
for (const name of names) {
  const text = readFileSync(join(src, `${name}.mmd`), 'utf8');
  const fontSize = Number((text.match(/^%%\s*fontSize:\s*(\d+)/m) || [])[1] || 14);
  const { svg, width, height } = render(text, { fontFamily: 'sans-serif', fontSize, padding: 8 });
  writeFileSync(join(out, `${name}.svg`), svg.replace('style="', `style="${theme}background:#fff;`));
  console.log(`${name}.svg  ${Math.ceil(width)} x ${Math.ceil(height)} px`);
}
