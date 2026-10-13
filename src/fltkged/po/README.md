# Translations

User-interface strings marked with `FLTKGED_TR` use the gettext domain
`fltkged`. Update the portable-object template with:

```sh
xgettext --language=C++ --keyword=FLTKGED_TR --from-code=UTF-8 \
  --output=fltkged.pot ../src/*.cpp
```

Compile a translation as `fltkged.mo` and install it as
`LOCALE/LC_MESSAGES/fltkged.mo` below the installation's `share/locale`
directory. For a relocatable or staged bundle, set `FLTKGED_LOCALE_PATH` to
the directory containing the locale subdirectories before starting `fltkged`.
Command names, database object names, Tcl, and terminal input are deliberately
not translated.
