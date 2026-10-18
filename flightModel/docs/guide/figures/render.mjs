// Renders figures/src/*.mmd (Mermaid) to figures/pdf/*.pdf for the guide.
//     node figures/render.mjs [name ...]
// Pele (npm install) makes the SVG; headless Chrome prints it to a vector PDF,
// because Pele colours the diagram with CSS variables.
import { render } from 'pele';
import { execFileSync } from 'node:child_process';
import { mkdirSync, readFileSync, readdirSync, writeFileSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { basename, dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';

const here = dirname(fileURLToPath(import.meta.url));
const src = join(here, 'src');
const out = join(here, 'pdf');
const chrome = process.env.CHROME || 'google-chrome';
const font = 'TeX Gyre Heros';

const theme = `
  --pele-bg: #ffffff; --pele-fg: #1a1a1a; --pele-muted: #5b6573;
  --pele-line: #1f3a5f; --pele-surface: #eaf2fb; --pele-surface-alt: #fff4e0;
  --pele-border: #2f6db5; --pele-accent: #c77700; --pele-font: '${font}', sans-serif;
  --pele-radius: 3px;`;

mkdirSync(out, { recursive: true });
const tmp = join(tmpdir(), `guide-figures-${process.pid}`);
mkdirSync(tmp, { recursive: true });

const names = process.argv.length > 2
  ? process.argv.slice(2)
  : readdirSync(src).filter((f) => f.endsWith('.mmd')).map((f) => basename(f, '.mmd'));

for (const name of names) {
  const text = readFileSync(join(src, `${name}.mmd`), 'utf8');
  const { svg, width, height } = render(text, { fontFamily: font, fontSize: 14, padding: 8 });
  const w = Math.ceil(width), h = Math.ceil(height);
  const html = `<!doctype html><meta charset="utf-8"><style>
    @page { size: ${w}px ${h}px; margin: 0 }
    html, body { margin: 0; padding: 0; background: #fff; overflow: hidden }
    svg { display: block }
    .pele { ${theme} }
  </style><div class="pele" style="width:${w}px;height:${h}px">${svg}</div>`;
  const page = join(tmp, `${name}.html`);
  writeFileSync(page, html);
  execFileSync(chrome, ['--headless=new', '--disable-gpu', '--no-pdf-header-footer',
    `--print-to-pdf=${join(out, `${name}.pdf`)}`, `file://${page}`], { stdio: 'ignore' });
  console.log(`${name}.pdf  ${w} x ${h} px`);
}
rmSync(tmp, { recursive: true, force: true });
