//! A desktop viewer for STEP files.
//!
//! ```text
//! cargo run --release -p golf_viewer [-- part.step]
//! ```
//!
//! Open a file with O (a file dialog), by dropping it on the window, or by
//! naming it on the command line. It's imported, healed and meshed in the
//! background, then shown with its colours; the window's title says what came
//! in. Drag with the left button to orbit, with the right to pan, scroll to
//! zoom, and press F to frame everything again.

mod load;

use std::path::PathBuf;
use std::sync::Mutex;
use std::sync::mpsc::Receiver;
use std::sync::mpsc::Sender;
use std::sync::mpsc::channel;

use bevy::prelude::*;
use bevy::window::FileDragAndDrop;
use golf_view::ShowScene;
use golf_view::ViewerPlugin;

use crate::load::Loaded;
use crate::load::load;

const TITLE: &str = "golf viewer";

fn main() -> AppExit {
    let (sender, receiver) = channel();
    let first = std::env::args().nth(1).map(PathBuf::from);
    if let Some(path) = &first {
        start(&sender, path.clone());
    }
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: match &first {
                    Some(path) => loading_title(path),
                    None => format!("{TITLE} (O to open a STEP file, or drop one here)"),
                },
                ..default()
            }),
            ..default()
        }))
        .add_plugins(ViewerPlugin)
        .insert_resource(Loads {
            sender,
            receiver: Mutex::new(receiver),
        })
        .add_systems(Update, (open_on_key, open_dropped, receive))
        .run()
}

/// Loads under way, reporting back over a channel.
#[derive(Resource)]
struct Loads {
    sender: Sender<Loaded>,
    receiver: Mutex<Receiver<Loaded>>,
}

/// Loads `path` on a thread of its own, so the window stays responsive.
fn start(sender: &Sender<Loaded>, path: PathBuf) {
    let sender = sender.clone();
    std::thread::spawn(move || {
        let _ = sender.send(load(&path));
    });
}

fn loading_title(path: &std::path::Path) -> String {
    format!("{TITLE}: loading {}…", file_name(path))
}

fn file_name(path: &std::path::Path) -> String {
    path.file_name().map_or_else(
        || path.display().to_string(),
        |n| n.to_string_lossy().into_owned(),
    )
}

fn set_title(windows: &mut Query<&mut Window>, title: String) {
    for mut window in windows {
        window.title.clone_from(&title);
    }
}

/// O asks for a file, on a thread: the dialog would block the window.
fn open_on_key(
    keys: Res<ButtonInput<KeyCode>>,
    loads: Res<Loads>,
    mut windows: Query<&mut Window>,
) {
    if !keys.just_pressed(KeyCode::KeyO) {
        return;
    }
    set_title(&mut windows, format!("{TITLE}: choosing a file…"));
    let sender = loads.sender.clone();
    std::thread::spawn(move || {
        let picked = rfd::FileDialog::new()
            .set_title("Open a STEP file")
            .add_filter("STEP", &["step", "stp", "STEP", "STP"])
            .pick_file();
        match picked {
            Some(path) => {
                let _ = sender.send(Loaded::Started(path.clone()));
                let _ = sender.send(load(&path));
            }
            None => {
                let _ = sender.send(Loaded::Cancelled);
            }
        }
    });
}

fn open_dropped(
    mut drops: MessageReader<FileDragAndDrop>,
    loads: Res<Loads>,
    mut windows: Query<&mut Window>,
) {
    for drop in drops.read() {
        if let FileDragAndDrop::DroppedFile { path_buf, .. } = drop {
            set_title(&mut windows, loading_title(path_buf));
            start(&loads.sender, path_buf.clone());
        }
    }
}

/// Shows what loads send back.
fn receive(loads: Res<Loads>, mut show: ResMut<ShowScene>, mut windows: Query<&mut Window>) {
    let Ok(receiver) = loads.receiver.lock() else {
        return;
    };
    while let Ok(loaded) = receiver.try_recv() {
        match loaded {
            Loaded::Started(path) => set_title(&mut windows, loading_title(&path)),
            Loaded::Cancelled => set_title(
                &mut windows,
                format!("{TITLE} (O to open a STEP file, or drop one here)"),
            ),
            Loaded::Failed(path, error) => set_title(
                &mut windows,
                format!("{TITLE}: {} failed: {error}", file_name(&path)),
            ),
            Loaded::Scene {
                path,
                scene,
                summary,
            } => {
                set_title(
                    &mut windows,
                    format!("{TITLE}: {} ({summary})", file_name(&path)),
                );
                show.0 = Some(scene);
            }
        }
    }
}
