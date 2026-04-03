/* ============================================================
   script.js — DNA Storage Encoder
   All logic and interactivity lives here.

   SECTIONS:
   1. Constants
   2. State variables
   3. Tab switching
   4. File upload handlers
   5. Encoding pipeline
   6. Decoding pipeline
   7. Display / render helpers
   8. Canvas helix drawing
   9. Utility functions
   ============================================================ */


/* ── 1. Constants ─────────────────────────────────────────── */

// Forward and reverse primers added to every oligo
const PRIMER_F = 'ATCGATCG';
const PRIMER_R = 'CGATCGAT';

// Max total bases per oligo (industry standard ~200)
const OLIGO_MAX = 200;

// Data bases per oligo = 200 - primer lengths
const OLIGO_DATA = OLIGO_MAX - PRIMER_F.length - PRIMER_R.length; // 184

// 2-bit → base mapping
const BIT_TO_BASE = {
  '00': 'A',
  '01': 'T',
  '10': 'G',
  '11': 'C'
};

// base → 2-bit mapping (reverse)
const BASE_TO_BIT = {
  'A': '00',
  'T': '01',
  'G': '10',
  'C': '11'
};

// Watson-Crick complement pairs (for the helix visual)
const COMPLEMENT = { A:'T', T:'A', G:'C', C:'G' };

// Colours for each base (used in canvas drawing)
const BASE_COLORS = {
  A: '#00c896',
  T: '#D85A30',
  G: '#7F77DD',
  C: '#378ADD'
};


/* ── 2. State variables ───────────────────────────────────── */

let currentTab    = 'text';  // which tab is open
let fullDNA       = '';      // raw DNA (no parity)
let fullDNAwithEC = '';      // DNA with parity bases inserted


/* ── 3. Tab switching ─────────────────────────────────────── */

function switchTab(tab) {
  currentTab = tab;

  const tabNames = ['text', 'file', 'decode'];

  tabNames.forEach((name, index) => {
    // show/hide tab content
    document.getElementById('tab-' + name).classList.toggle('hidden', name !== tab);

    // highlight active tab button
    document.querySelectorAll('.tab')[index].classList.toggle('active', name === tab);
  });
}


/* ── 4. File upload handlers ──────────────────────────────── */

// When user picks a file via the input button
document.getElementById('file-input').addEventListener('change', function(e) {
  const file = e.target.files[0];
  if (file) {
    document.getElementById('upload-label').textContent =
      '✓ ' + file.name + '  (' + formatBytes(file.size) + ')';
  }
});

// Drag over: highlight the zone
const dropZone = document.getElementById('drop-zone');
dropZone.addEventListener('dragover', function(e) {
  e.preventDefault();
  dropZone.style.background = 'var(--green-dim)';
});

// Drag leave: un-highlight
dropZone.addEventListener('dragleave', function() {
  dropZone.style.background = '';
});

// Drop: load the dropped file
dropZone.addEventListener('drop', function(e) {
  e.preventDefault();
  dropZone.style.background = '';
  const file = e.dataTransfer.files[0];
  if (file) {
    document.getElementById('file-input').files = e.dataTransfer.files;
    document.getElementById('upload-label').textContent =
      '✓ ' + file.name + '  (' + formatBytes(file.size) + ')';
  }
});


/* ── 5. Encoding pipeline ─────────────────────────────────── */

/*
  Main encode() — called by both the text and file buttons.
  Reads input, converts to bytes, then calls runEncode().
*/
function encode() {

  if (currentTab === 'text') {
    const text = document.getElementById('text-input').value;
    if (!text.trim()) return;

    // TextEncoder converts a JS string → Uint8Array of UTF-8 bytes
    const bytes = new TextEncoder().encode(text);
    runEncode(bytes, text.length + ' chars');

  } else {
    // File tab
    const file = document.getElementById('file-input').files[0];
    if (!file) return;

    // Limit to 200 KB so the browser doesn't freeze on huge files
    const MAX_BYTES = 200 * 1024;
    const reader = new FileReader();

    reader.onload = function(e) {
      let bytes = new Uint8Array(e.target.result);
      const wasTruncated = bytes.length > MAX_BYTES;
      if (wasTruncated) bytes = bytes.slice(0, MAX_BYTES);

      const label = (wasTruncated ? 'Preview 200 KB of ' : '') + file.name;
      runEncode(bytes, label);
    };

    reader.readAsArrayBuffer(file);
  }
}

/*
  runEncode(bytes, label)
  Full pipeline:
    bytes → binary string → DNA → add parity → chunk into oligos → display
*/
function runEncode(bytes, label) {

  // Step 1: bytes → binary string  e.g. "01001000 01101001 …"
  const binaryString = bytesToBinary(bytes);

  // Step 2: binary → raw DNA  e.g. "TAGATGT…"
  const dnaRaw = binaryToDNA(binaryString);

  // Step 3: insert parity bases for error detection
  const { sequence: dnaWithParity, parityPositions } = addParityBases(dnaRaw);

  // Step 4: split into oligos (each has primers + max 184 data bases)
  const oligos = splitIntoOligos(dnaWithParity);

  // Save to state so other functions can use them
  fullDNA       = dnaRaw;
  fullDNAwithEC = dnaWithParity;

  // Step 5: show everything on screen
  displayResults(bytes, binaryString, dnaRaw, dnaWithParity, parityPositions, oligos, label);
}

/* -- Encoding helper functions -- */

// Convert bytes to a long binary string
// e.g. [72] → "01001000"
function bytesToBinary(bytes) {
  let binary = '';
  for (let i = 0; i < bytes.length; i++) {
    binary += bytes[i].toString(2).padStart(8, '0');
  }
  return binary;
}

// Convert binary string to DNA bases using 2-bit mapping
// e.g. "01001000" → "TAGA"
function binaryToDNA(binary) {
  // Pad to even length if needed
  if (binary.length % 2 !== 0) binary += '0';

  let dna = '';
  for (let i = 0; i < binary.length; i += 2) {
    const bits = binary[i] + binary[i + 1];
    dna += BIT_TO_BASE[bits];
  }
  return dna;
}

/*
  addParityBases(dna)
  After every 8 data bases (= 1 byte), insert 1 parity base.
  Parity rule (even parity on G/C count):
    - Count how many G or C bases are in the 8-base block
    - If count is EVEN  → parity base = A
    - If count is ODD   → parity base = T
  Returns { sequence, parityPositions }
*/
function addParityBases(dna) {
  let result = '';
  let position = 0;
  const parityPositions = [];

  for (let i = 0; i < dna.length; i += 8) {
    const block = dna.slice(i, i + 8);
    result += block;
    position += block.length;

    // Count G and C in block
    const gcCount = [...block].filter(b => b === 'G' || b === 'C').length;
    const parityBase = (gcCount % 2 === 0) ? 'A' : 'T';

    result += parityBase;
    parityPositions.push(position); // record where the parity base is
    position += 1;
  }

  return { sequence: result, parityPositions };
}

/*
  splitIntoOligos(dna)
  Split the full sequence into chunks of OLIGO_DATA bases,
  then wrap each chunk with PRIMER_F and PRIMER_R.
*/
function splitIntoOligos(dna) {
  const oligos = [];

  for (let i = 0; i < dna.length; i += OLIGO_DATA) {
    const chunk    = dna.slice(i, i + OLIGO_DATA);
    const full     = PRIMER_F + chunk + PRIMER_R;
    oligos.push({
      index:    oligos.length + 1,
      data:     chunk,
      full:     full,
      length:   full.length
    });
  }

  return oligos;
}


/* ── 6. Decoding pipeline ─────────────────────────────────── */

/*
  decodeDNA()
  Called by the Decode button.
  Reads DNA input, optionally strips parity, converts back to text.
*/
function decodeDNA() {
  // Clean input: uppercase and keep only valid bases
  let raw = document.getElementById('dna-decode-input').value
    .trim().toUpperCase().replace(/[^ATGC]/g, '');

  if (!raw) return;

  const hasParity = document.getElementById('has-parity').checked;
  let dataDNA = '';
  let parityErrors = [];

  if (hasParity) {
    // Strip parity bases and check them
    for (let i = 0; i < raw.length; i += 9) {
      const block      = raw.slice(i, i + 8);
      const parityBase = raw[i + 8];

      // Add data block
      dataDNA += block;

      if (!parityBase) break; // last block might be incomplete

      // Verify parity
      const gcCount    = [...block].filter(b => b === 'G' || b === 'C').length;
      const expectedPB = (gcCount % 2 === 0) ? 'A' : 'T';
      if (parityBase !== expectedPB) {
        parityErrors.push(Math.floor(i / 9) + 1); // record block number
      }
    }
  } else {
    dataDNA = raw;
  }

  // Convert DNA → binary → bytes
  const bits  = [...dataDNA].map(b => BASE_TO_BIT[b] || '00').join('');
  const bytes = [];
  for (let i = 0; i + 7 < bits.length; i += 8) {
    bytes.push(parseInt(bits.slice(i, i + 8), 2));
  }

  const uint8 = new Uint8Array(bytes);

  // Try to decode as UTF-8 text
  let outputText = '';
  try {
    outputText = new TextDecoder('utf-8', { fatal: true }).decode(uint8);
  } catch(e) {
    // Not valid text — show hex preview
    const hex = [...uint8.slice(0, 64)].map(b => b.toString(16).padStart(2,'0')).join(' ');
    outputText = '[Binary file — not displayable as text]\n\nFirst 64 bytes (hex):\n' + hex;
  }

  // Show results
  document.getElementById('decode-out').textContent = outputText;

  const statusEl = document.getElementById('decode-status');
  const totalBlocks = Math.floor(raw.length / 9);

  if (hasParity && parityErrors.length === 0) {
    statusEl.innerHTML =
      '<span style="color:var(--green);font-size:12px;">✓ Parity check passed — all ' + totalBlocks + ' blocks OK</span>';
  } else if (hasParity && parityErrors.length > 0) {
    statusEl.innerHTML =
      '<span style="color:#D85A30;font-size:12px;">✗ Parity errors in ' + parityErrors.length +
      ' block(s): ' + parityErrors.slice(0, 8).join(', ') + '</span>';
  } else {
    statusEl.innerHTML =
      '<span style="color:var(--text-dim);font-size:12px;">ℹ Decoded without parity check · ' +
      dataDNA.length + ' bases → ' + bytes.length + ' bytes</span>';
  }

  document.getElementById('decode-results').style.display = 'block';
}


/* ── 7. Display / render helpers ──────────────────────────── */

/*
  displayResults(...)
  Populates every section of the results panel.
*/
function displayResults(bytes, binary, dnaRaw, dnaEC, parityPositions, oligos, label) {

  // Success banner
  document.getElementById('success-msg').textContent =
    label + ' encoded → ' + formatNumber(dnaEC.length) + ' DNA bases';

  // Stat cards
  document.getElementById('s-bytes').textContent    = formatNumber(bytes.length);
  document.getElementById('s-bases').textContent    = formatNumber(dnaRaw.length);
  document.getElementById('s-bases-ec').textContent = formatNumber(dnaEC.length);
  document.getElementById('s-oligos').textContent   = formatNumber(oligos.length);

  // GC content
  renderGCBar(dnaRaw);

  // Binary preview (first 512 bits)
  document.getElementById('binary-out').textContent =
    binary.slice(0, 512) + (binary.length > 512 ? '…' : '');

  // DNA sequence preview (first 300 bases, with parity highlighted)
  renderDNAPreview(dnaEC, parityPositions);

  // Error correction explanation
  renderECExplain(dnaRaw);

  // Oligo cards
  renderOligos(oligos);

  // Show results panel
  document.getElementById('results').style.display = 'block';

  // Draw helix canvas
  drawHelix(dnaRaw.slice(0, 200));
}

// Render the GC% bar and base breakdown
function renderGCBar(dna) {
  const counts = { A:0, T:0, G:0, C:0 };
  for (const base of dna) counts[base]++;

  const total  = dna.length || 1;
  const gcPct  = Math.round((counts.G + counts.C) / total * 100);
  const gcColor = (gcPct > 70 || gcPct < 30) ? '#D85A30' : '#00c896';

  document.getElementById('gc-bar').style.width      = gcPct + '%';
  document.getElementById('gc-bar').style.background = gcColor;
  document.getElementById('gc-pct-lbl').textContent  = 'GC: ' + gcPct + '%';

  document.getElementById('at-breakdown').textContent =
    'A: ' + Math.round(counts.A / total * 100) + '%  ' +
    'T: ' + Math.round(counts.T / total * 100) + '%';

  document.getElementById('gc-breakdown').textContent =
    'G: ' + Math.round(counts.G / total * 100) + '%  ' +
    'C: ' + Math.round(counts.C / total * 100) + '%';

  // GC warning banner
  const warnEl  = document.getElementById('gc-warn');
  const warnMsg = document.getElementById('gc-warn-msg');

  if (gcPct > 70) {
    warnMsg.textContent = 'High GC (' + gcPct + '%) — above 70% makes synthesis harder and risks secondary structures.';
    warnEl.classList.remove('hidden');
  } else if (gcPct < 30) {
    warnMsg.textContent = 'Low GC (' + gcPct + '%) — below 30% reduces stability and increases synthesis error rate.';
    warnEl.classList.remove('hidden');
  } else {
    warnEl.classList.add('hidden');
  }
}

// Render coloured DNA preview with parity bases in orange
function renderDNAPreview(dnaEC, parityPositions) {
  const paritySet = new Set(parityPositions);
  const preview   = dnaEC.slice(0, 300);
  let html = '';

  for (let i = 0; i < preview.length; i++) {
    const base = preview[i];
    if (paritySet.has(i)) {
      html += '<span class="parity">' + base + '</span>';
    } else {
      html += '<span class="' + base + '">' + base + '</span>';
    }
  }

  if (dnaEC.length > 300) {
    html += '<span style="color:var(--text-muted)">…+' + formatNumber(dnaEC.length - 300) + ' more</span>';
  }

  document.getElementById('dna-out').innerHTML = html;
}

// Show the first few parity blocks explained
function renderECExplain(dnaRaw) {
  const lines = [];
  for (let i = 0; i < Math.min(dnaRaw.length, 32); i += 8) {
    const block   = dnaRaw.slice(i, i + 8);
    const gcCount = [...block].filter(b => b === 'G' || b === 'C').length;
    const parity  = (gcCount % 2 === 0) ? 'A' : 'T';
    lines.push(
      'Block ' + (i / 8 + 1) + ': [' + block + ']  GC=' + gcCount +
      '  parity=' + parity + '  (' + (gcCount % 2 === 0 ? 'even→A' : 'odd→T') + ')'
    );
  }
  document.getElementById('ec-explain').textContent = lines.join('\n');
}

// Render oligo cards
function renderOligos(oligos) {
  const container = document.getElementById('oligo-wrap');
  container.innerHTML = '';

  const showCount = Math.min(oligos.length, 20);

  for (let o = 0; o < showCount; o++) {
    const oligo = oligos[o];
    const div   = document.createElement('div');
    div.className = 'oligo';

    // Build coloured sequence HTML
    let seqHTML = '<span class="primer">' + PRIMER_F + '</span>';

    const dataPreview = oligo.data.slice(0, 55);
    for (const base of dataPreview) {
      seqHTML += '<span class="' + base + '">' + base + '</span>';
    }
    if (oligo.data.length > 55) {
      seqHTML += '<span style="color:var(--text-muted)">…+' + (oligo.data.length - 55) + ' more</span>';
    }
    seqHTML += '<span class="primer">' + PRIMER_R + '</span>';

    div.innerHTML =
      '<div class="oligo-head">' +
        '<span class="oligo-id">oligo #' + oligo.index + '</span>' +
        '<span class="oligo-id">' + oligo.length + ' bases · ' + oligo.data.length + ' data</span>' +
      '</div>' +
      '<div class="oligo-seq">' + seqHTML + '</div>';

    container.appendChild(div);
  }

  // "and N more" footer
  if (oligos.length > 20) {
    const more = document.createElement('div');
    more.style.cssText = 'font-size:10px;color:var(--text-muted);padding:5px 0;';
    more.textContent   = '… and ' + (oligos.length - 20) + ' more oligos';
    container.appendChild(more);
  }
}


/* ── 8. Canvas helix drawing ──────────────────────────────── */

/*
  drawHelix(sequence)
  Draws a sine-wave double helix on the canvas.
  Top strand = sense strand (bases from sequence)
  Bottom strand = antisense strand (complementary bases)
  Vertical connectors = hydrogen bonds between base pairs
*/
function drawHelix(sequence) {
  const canvas  = document.getElementById('dna-canvas');
  const W       = canvas.offsetWidth || 800;
  const H       = 130;

  // Scale for high-DPI / retina screens
  canvas.width  = W * devicePixelRatio;
  canvas.height = H * devicePixelRatio;
  canvas.style.height = H + 'px';

  const ctx = canvas.getContext('2d');
  ctx.scale(devicePixelRatio, devicePixelRatio);
  ctx.clearRect(0, 0, W, H);

  const n         = Math.min(sequence.length, Math.floor(W / 5)); // bases that fit
  const baseWidth = W / n;
  const centerY   = H / 2;
  const amplitude = 26; // wave height
  const period    = 20; // bases per full wave cycle

  // Draw backbone strands (connecting lines between adjacent bases)
  for (let i = 0; i < n - 1; i++) {
    const x1 = i * baseWidth + baseWidth / 2;
    const x2 = (i + 1) * baseWidth + baseWidth / 2;

    // y positions on the sine wave
    const y1_top = centerY - amplitude * Math.sin(i * Math.PI * 2 / period);
    const y2_top = centerY - amplitude * Math.sin((i + 1) * Math.PI * 2 / period);
    const y1_bot = centerY + amplitude * Math.sin(i * Math.PI * 2 / period);
    const y2_bot = centerY + amplitude * Math.sin((i + 1) * Math.PI * 2 / period);

    // Top backbone
    ctx.beginPath();
    ctx.strokeStyle = BASE_COLORS[sequence[i]] + '88';
    ctx.lineWidth   = 1.5;
    ctx.moveTo(x1, y1_top);
    ctx.lineTo(x2, y2_top);
    ctx.stroke();

    // Bottom backbone (complementary strand)
    ctx.beginPath();
    ctx.strokeStyle = BASE_COLORS[COMPLEMENT[sequence[i]]] + '88';
    ctx.lineWidth   = 1.5;
    ctx.moveTo(x1, y1_bot);
    ctx.lineTo(x2, y2_bot);
    ctx.stroke();
  }

  // Draw base pair connectors and dots
  for (let i = 0; i < n; i++) {
    const x   = i * baseWidth + baseWidth / 2;
    const y1  = centerY - amplitude * Math.sin(i * Math.PI * 2 / period);
    const y2  = centerY + amplitude * Math.sin(i * Math.PI * 2 / period);

    // Hydrogen bond connector (vertical line between the two strands)
    ctx.beginPath();
    ctx.strokeStyle = BASE_COLORS[sequence[i]] + '33';
    ctx.lineWidth   = 0.8;
    ctx.moveTo(x, y1);
    ctx.lineTo(x, y2);
    ctx.stroke();

    // Top strand dot (sense base)
    ctx.beginPath();
    ctx.fillStyle = BASE_COLORS[sequence[i]];
    ctx.arc(x, y1, 3, 0, Math.PI * 2);
    ctx.fill();

    // Bottom strand dot (complementary base)
    ctx.beginPath();
    ctx.fillStyle = BASE_COLORS[COMPLEMENT[sequence[i]]];
    ctx.arc(x, y2, 3, 0, Math.PI * 2);
    ctx.fill();
  }
}


/* ── 9. Utility functions ─────────────────────────────────── */

// Copy full DNA sequence (with parity) to clipboard
function copyDNA() {
  navigator.clipboard.writeText(fullDNAwithEC).then(function() {
    const btn = document.getElementById('copy-btn');
    btn.textContent = 'copied!';
    setTimeout(function() { btn.textContent = 'copy full'; }, 1500);
  });
}

// Format a number with commas  e.g. 1000000 → "1,000,000"
function formatNumber(n) {
  return n.toLocaleString();
}

// Format bytes into readable size  e.g. 2048 → "2.0 KB"
function formatBytes(bytes) {
  if (bytes < 1024)        return bytes + ' B';
  if (bytes < 1024 * 1024) return (bytes / 1024).toFixed(1) + ' KB';
  return (bytes / (1024 * 1024)).toFixed(2) + ' MB';
}
